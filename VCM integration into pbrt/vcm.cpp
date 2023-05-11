
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// integrators/vcm.cpp*
#include "integrators/vcm.h"
#include "film.h"
#include "filters/box.h"
#include "integrator.h"
#include "lightdistrib.h"
#include "paramset.h"
#include "progressreporter.h"
#include "sampler.h"
#include "stats.h"
#include "error.h"

namespace pbrt {

STAT_PERCENT("Integrator/Zero-radiance paths", zeroRadiancePaths, totalPaths);
STAT_INT_DISTRIBUTION("Integrator/Path length", pathLength);

// VCM Forward Declarations
int VCMRandomWalk(const Scene &scene, RayDifferential ray, Sampler &sampler,
                  MemoryArena &arena, VCMConstant vcmConstant, Spectrum beta,
                  Float pdf, int maxDepth, TransportMode mode, VCMVertex *path);

// VCM Grid Utility Functions
// [directly adopted from pbrt]
static bool ToGrid(const Point3f &p, const Bounds3f &bounds,
                   const int gridRes[3], Point3i *pi) {
    bool inBounds = true;
    Vector3f pg = bounds.Offset(p);
    for (int i = 0; i < 3; ++i) {
        (*pi)[i] = (int)(gridRes[i] * pg[i]);
        inBounds &= ((*pi)[i] >= 0 && (*pi)[i] < gridRes[i]);
        (*pi)[i] = Clamp((*pi)[i], 0, gridRes[i] - 1);
    }
    return inBounds;
}

// [directly adopted from pbrt]
inline unsigned int hash(const Point3i &p, int hashSize) {
    return (unsigned int)((p.x * 73856093) ^ (p.y * 19349663) ^
                          (p.z * 83492791)) %
           hashSize;
}

// VCM Utility Functions
// [directly adopted from pbrt]
Float VCMCorrectShadingNormal(const SurfaceInteraction &isect, const Vector3f &wo,
                           const Vector3f &wi, TransportMode mode) {
    if (mode == TransportMode::Importance) {
        Float num = AbsDot(wo, isect.shading.n) * AbsDot(wi, isect.n);
        Float denom = AbsDot(wo, isect.n) * AbsDot(wi, isect.shading.n);
        // wi is occasionally perpendicular to isect.shading.n; this is
        // fine, but we don't want to return an infinite or NaN value in
        // that case.
        if (denom == 0) return 0;
        return num / denom;
    } else
        return 1;
}

// [adopted from pbrt and modified]
int VCMGenerateCameraSubpath(const Scene &scene, Sampler &sampler,
                             MemoryArena &arena, VCMConstant vcmConstant,
                             int maxDepth, const Camera &camera,
                             const Point2f &pFilm, VCMVertex *path) {
    if (maxDepth == 0) return 0;

    ProfilePhase _(Prof::BDPTGenerateSubpath);

    // Sample initial ray for camera subpath
    CameraSample cameraSample;
    cameraSample.pFilm = pFilm;
    cameraSample.time = sampler.Get1D();
    cameraSample.pLens = sampler.Get2D();

    // Generate ray
    RayDifferential ray;
    Spectrum beta = camera.GenerateRayDifferential(cameraSample, &ray);
    ray.ScaleDifferentials(1 / std::sqrt(sampler.samplesPerPixel));

    // Generate first vertex on camera subpath (on camera)
    Float pdfPos, pdfDir;
    camera.Pdf_We(ray, &pdfPos, &pdfDir);
    path[0] = VCMVertex::CreateCamera(&camera, ray, beta, pdfPos, pdfDir);

    //VLOG(2) << "Starting camera subpath. Ray: " << ray << ", beta " << beta
    //        << ", pdfPos " << pdfPos << ", pdfDir " << pdfDir;

    // Start random walk
    int nVertices =
        VCMRandomWalk(scene, ray, sampler, arena, vcmConstant, beta, pdfDir,
                      maxDepth - 1, TransportMode::Radiance, path + 1);
    return nVertices + 1;
}

// [adopted from pbrt and modified]
int VCMGenerateLightSubpath(
    const Scene &scene, Sampler &sampler, MemoryArena &arena,
    VCMConstant vcmConstant, int maxDepth, Float time,
    const Distribution1D &lightDistr,
    const std::unordered_map<const Light *, size_t> &lightToIndex,
    VCMVertex *path) {
    if (maxDepth == 0) return 0;

    ProfilePhase _(Prof::BDPTGenerateSubpath);

    // Sample a light source
    Float lightPdf;
    int lightNum = lightDistr.SampleDiscrete(sampler.Get1D(), &lightPdf);
    const std::shared_ptr<Light> &light = scene.lights[lightNum];

    // Sample initial ray for light subpath
    RayDifferential ray;
    Normal3f nLight;
    Float pdfPos, pdfDir;
    Spectrum Le = light->Sample_Le(sampler.Get2D(), sampler.Get2D(), time, &ray,
                                   &nLight, &pdfPos, &pdfDir);
    if (pdfPos == 0 || pdfDir == 0 || Le.IsBlack()) return 0;

    // Generate first vertex on light subpath (on light source)
    path[0] =
        VCMVertex::CreateLight(light.get(), ray, nLight, Le, pdfPos * lightPdf,
                               lightPdf, pdfPos, pdfDir);
    Spectrum beta = Le * AbsDot(nLight, ray.d) / (lightPdf * pdfPos * pdfDir);

    //VLOG(2) << "Starting light subpath. Ray: " << ray << ", Le " << Le <<
    //    ", beta " << beta << ", pdfPos " << pdfPos << ", pdfDir " << pdfDir;

    // Start random walk
    int nVertices =
        VCMRandomWalk(scene, ray, sampler, arena, vcmConstant, beta, pdfDir,
                      maxDepth - 1, TransportMode::Importance, path + 1);

    // Correct subpath sampling densities for infinite area lights
    if (path[0].IsInfiniteLight()) {
        // Set spatial density of _path[1]_ for infinite area light
        if (nVertices > 0) {
            path[1].pdfFwd = pdfPos;
            if (path[1].IsOnSurface())
                path[1].pdfFwd *= AbsDot(ray.d, path[1].ng());
        }

        // Set spatial density of _path[0]_ for infinite area light
        path[0].pdfFwd =
            VCMInfiniteLightDensity(scene, lightDistr, lightToIndex, ray.d);
    }
    return nVertices + 1;
}

// [adopted from pbrt and modified]
int VCMRandomWalk(const Scene &scene, RayDifferential ray, Sampler &sampler,
                  MemoryArena &arena, VCMConstant vcmConstant, Spectrum beta,
                  Float pdf, int maxDepth, TransportMode mode,
                  VCMVertex *path) {
    if (maxDepth == 0) return 0;

    // This bounces variable corresponds to the number of edges in the path
    // i.e., v0 -> v1 gives bounces = 1
    int bounces = 0;

    // Declare variables for forward and reverse probability densities
    Float pdfFwd = pdf, pdfRev = 0;

    while (true) {
        // Attempt to create the next subpath vertex in _path_

        //VLOG(2) << "Random walk. Bounces " << bounces << ", beta " << beta <<
        //    ", pdfFwd " << pdfFwd << ", pdfRev " << pdfRev;

        // Trace a ray
        SurfaceInteraction isect;
        bool foundIntersection = scene.Intersect(ray, &isect);

        // Sample the medium, if any
        MediumInteraction mi;
        if (ray.medium) beta *= ray.medium->Sample(ray, sampler, arena, &mi);
        if (beta.IsBlack()) break;

        VCMVertex &vertex = path[bounces], &prev = path[bounces - 1];
        VCMVertex *prev2 = (bounces - 2 >= -1) ? &path[bounces - 2] : nullptr;

        if (mi.IsValid()) {
            // Record medium interaction in _path_ and compute forward density
            vertex =
                VCMVertex::CreateMedium(mi, beta, pdfFwd, prev);
            if (++bounces >= maxDepth) break;

            // Sample direction and compute reverse density at preceding vertex
            Vector3f wi;
            pdfFwd = pdfRev = mi.phase->Sample_p(-ray.d, &wi, sampler.Get2D());
            ray = mi.SpawnRay(wi);

        } else {
            // Handle surface interaction for path generation

            // Capture escaped rays when tracing from the camera
            if (!foundIntersection) {
                if (mode == TransportMode::Radiance) {
                    vertex = VCMVertex::CreateLight(VCMEndpointInteraction(ray),
                                                    beta, pdfFwd);
                    ++bounces;
                }
                break;
            }

            // Compute scattering functions for _mode_ and skip over medium
            // boundaries
            isect.ComputeScatteringFunctions(ray, arena, true, mode);
            if (!isect.bsdf) {
                ray = isect.SpawnRay(ray.d);
                continue;
            }

            // Initialize _vertex_ with surface intersection information
            vertex = VCMVertex::CreateSurface(isect, beta, pdfFwd, prev);
            if (++bounces >= maxDepth) break;

            // Sample BSDF at current vertex and get forward pdf for next vertex
            Vector3f wi, wo = isect.wo;
            BxDFType type;
            Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdfFwd,
                                              BSDF_ALL, &type);
            if (f.IsBlack() || pdfFwd == 0.f) break;
            beta *= f * AbsDot(wi, isect.shading.n) / pdfFwd;

            //VLOG(2) << "Random walk sampled dir " << wi << " f: " << f <<
            //    ", pdfFwd: " << pdfFwd;
            //VLOG(2) << "Random walk beta now " << beta;

            // Compute reverse probability
            pdfRev = isect.bsdf->Pdf(wi, wo, BSDF_ALL);

            // Set pdfs to zero for specular reflections
            if (type & BSDF_SPECULAR) {
                vertex.delta = true;
                pdfRev = pdfFwd = 0;
            }

            // Account for shading normal
            beta *= VCMCorrectShadingNormal(isect, wo, wi, mode);

            //VLOG(2) << "Random walk beta after shading normal correction " << beta;

            // Spawn a new way to continue random walk
            ray = isect.SpawnRay(wi);
        }

        // Compute reverse area density at preceding vertex
        prev.pdfRev = vertex.ConvertDensity(pdfRev, prev);
    }
    return bounces;
}

// [directly adopted from pbrt]
Spectrum VCMG(const Scene &scene, Sampler &sampler, const VCMVertex &v0,
           const VCMVertex &v1) {
    Vector3f d = v0.p() - v1.p();
    Float g = 1 / d.LengthSquared();
    d *= std::sqrt(g);
    if (v0.IsOnSurface()) g *= AbsDot(v0.ns(), d);
    if (v1.IsOnSurface()) g *= AbsDot(v1.ns(), d);
    VisibilityTester vis(v0.GetInteraction(), v1.GetInteraction());
    return g * vis.Tr(scene, sampler);
}

// [adopted from pbrt and substantially modified]
Float VCMMISWeight(
    const Scene &scene, VCMConstant vcmConstant, VCMPathStrategy pathStrategy,
    VCMVertex *lightVertices, VCMVertex *cameraVertices, VCMVertex &sampled,
    int s, int t, const Distribution1D &lightPdf,
    const std::unordered_map<const Light *, size_t> &lightToIndex) {
    if (s + t == 2) return 1;
    Float sumRi = 0;
    // Define helper function _remap0_ that deals with Dirac delta functions
    auto remap0 = [](Float f) -> Float { return f != 0 ? f : 1; };

    // Temporarily update vertex properties for current strategy

    // Look up connection vertices and their predecessors
    VCMVertex *qs = s > 0 ? &lightVertices[s - 1] : nullptr,
           *pt = t > 0 ? &cameraVertices[t - 1] : nullptr,
           *qsMinus = s > 1 ? &lightVertices[s - 2] : nullptr,
           *ptMinus = t > 1 ? &cameraVertices[t - 2] : nullptr;

    // Update sampled vertex for $s=1$ or $t=1$ strategy
    VCMScopedAssignment<VCMVertex> a1;
    if (s == 1)
        a1 = {qs, sampled};
    else if (t == 1)
        a1 = {pt, sampled};

    // Mark connection vertices as non-degenerate
    VCMScopedAssignment<bool> a2, a3;
    if (pt) a2 = {&pt->delta, false};
    if (qs) a3 = {&qs->delta, false};

    // Update reverse density of vertex $\pt{}_{t-1}$
    VCMScopedAssignment<Float> a4;
    if (pt)
        a4 = {&pt->pdfRev, s > 0 ? qs->Pdf(scene, qsMinus, *pt)
                                 : pt->PdfLightOrigin(scene, *ptMinus, lightPdf,
                                                      lightToIndex)};

    // Update reverse density of vertex $\pt{}_{t-2}$
    VCMScopedAssignment<Float> a5;
    if (ptMinus)
        a5 = {&ptMinus->pdfRev, s > 0 ? pt->Pdf(scene, qs, *ptMinus)
                                      : pt->PdfLight(scene, *ptMinus)};

    // Update reverse density of vertices $\pq{}_{s-1}$ and $\pq{}_{s-2}$
    VCMScopedAssignment<Float> a6;
    if (qs) a6 = {&qs->pdfRev, pt->Pdf(scene, ptMinus, *qs)};
    VCMScopedAssignment<Float> a7;
    if (qsMinus) a7 = {&qsMinus->pdfRev, qs->Pdf(scene, pt, *qsMinus)};

    // Weight computation
    Float wLight = 0;
    Float wEye = 0;
    if (pathStrategy == VC) {
        // Vertex connection

        Float wVCLightVCTermSum = 0;
        Float wVCLightVCTermProduct = 1;
        for (int j = s - 1; j >= 0; j--) {
            wVCLightVCTermProduct *= remap0(lightVertices[j].pdfRev) /
                                     remap0(lightVertices[j].pdfFwd);

            bool deltaLightVertex = j > 0 ? lightVertices[j - 1].delta
                                          : lightVertices[0].IsDeltaLight();

            if (!lightVertices[j].delta && !deltaLightVertex) {
                wVCLightVCTermSum += wVCLightVCTermProduct;
            }
        }

        Float wVCLightVMTermSum = 0;
        Float wVCLightVMTermProduct = 1;
        for (int j = s; j >= 2; j--) {
            if (j != s) {
                wVCLightVMTermProduct *= remap0(lightVertices[j].pdfRev) /
                                         remap0(lightVertices[j].pdfFwd);
            }

            bool deltaLightVertex =
                j < s ? lightVertices[j].delta : lightVertices[j - 1].delta;

            if (!deltaLightVertex && !lightVertices[j - 1].delta) {
                wVCLightVMTermSum +=
                    remap0(lightVertices[j - 1].pdfRev) * wVCLightVMTermProduct;
            }
        }

        Float wVCEyeVCTermSum = 0;
        Float wVCEyeVCTermProduct = 1;
        for (int j = t - 1; j > 0; j--) {
            wVCEyeVCTermProduct *= remap0(cameraVertices[j].pdfRev) /
                                   remap0(cameraVertices[j].pdfFwd);

            if (!cameraVertices[j].delta && !cameraVertices[j - 1].delta) {
                wVCEyeVCTermSum += wVCEyeVCTermProduct;
            }
        }

        Float wVCEyeVMTermSum = 0;
        Float wVCEyeVMTermProduct = 1;
        for (int j = t; j >= 2; j--) {
            if (j != t) {
                wVCEyeVMTermProduct *= remap0(cameraVertices[j].pdfRev) /
                                       remap0(cameraVertices[j].pdfFwd);
            }

            bool deltaCameraVertex =
                j < t ? cameraVertices[j].delta : cameraVertices[j - 1].delta;

            if (!deltaCameraVertex && !cameraVertices[j - 1].delta) {
                wVCEyeVMTermSum +=
                    remap0(cameraVertices[j - 1].pdfRev) * wVCEyeVMTermProduct;
            }
        }
        
        wLight = wVCLightVCTermSum + vcmConstant.etaVCM * wVCLightVMTermSum;
        wEye = wVCEyeVCTermSum + vcmConstant.etaVCM * wVCEyeVMTermSum;
		
    } else {
        // Vertex merging
        
        Float wVMLightVCTermSum = 0;
        Float wVMLightVCTermProduct = 1;
        for (int j = s - 1; j >= 0; j--) {
            if (j != s - 1) {
                wVMLightVCTermProduct *= remap0(lightVertices[j].pdfRev) /
                                         remap0(lightVertices[j].pdfFwd);
            }

            bool deltaLightVertex = j > 0 ? lightVertices[j - 1].delta
                                          : lightVertices[0].IsDeltaLight();

            if (!lightVertices[j].delta && !deltaLightVertex) {
                wVMLightVCTermSum += wVMLightVCTermProduct;
            }
        }

        Float wVMLightVMTermSum = 0;
        Float wVMLightVMTermProduct = 1;
        for (int j = s - 1; j >= 2; j--) {
            wVMLightVMTermProduct *= remap0(lightVertices[j - 1].pdfRev) /
                                     remap0(lightVertices[j].pdfFwd);

            if (!lightVertices[j].delta && !lightVertices[j - 1].delta) {
                wVMLightVMTermSum += wVMLightVMTermProduct;
            }
        }

        Float wVMEyeVCTermSum = 0;
        Float wVMEyeVCTermProduct = 1;
        for (int j = t - 1; j > 0; j--) {
            if (t != t - 1) {
                wVMEyeVCTermProduct *= remap0(cameraVertices[j].pdfRev) /
                                       remap0(cameraVertices[j].pdfFwd);
            }

            if (!cameraVertices[j].delta && !cameraVertices[j - 1].delta) {
                wVMEyeVCTermSum += wVMEyeVCTermProduct;
            }
        }

        Float wVMEyeVMTermSum = 0;
        Float wVMEyeVMTermProduct = 1;
        for (int j = t - 1; j >= 2; j--) {
            wVMEyeVMTermProduct *= remap0(cameraVertices[j - 1].pdfRev) /
                                   remap0(cameraVertices[j].pdfFwd);

            if (!cameraVertices[j].delta && !cameraVertices[j - 1].delta) {
                wVMEyeVMTermSum += wVMEyeVMTermProduct;
            }
        }

        wLight = vcmConstant.invEtaVCM * 1.f / lightVertices[s - 1].pdfFwd *
                     wVMLightVCTermSum +
                 wVMLightVMTermSum;
        wEye = vcmConstant.invEtaVCM * 1.f / cameraVertices[t - 1].pdfFwd *
                   wVMEyeVCTermSum +
               wVMEyeVMTermSum;
    }
	
    return 1.f / (wLight + 1 + wEye);
}

// VCM Method Definitions
// [Adopted from pbrt but rewritten completely]
void VCMIntegrator::Render(const Scene &scene) {
    std::unique_ptr<LightDistribution> lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);

    // Compute a reverse mapping from light pointers to offsets into the
    // scene lights vector (and, equivalently, offsets into
    // lightDistr). Added after book text was finalized; this is critical
    // to reasonable performance with 100s+ of light sources.
    std::unordered_map<const Light *, size_t> lightToIndex;
    for (size_t i = 0; i < scene.lights.size(); ++i)
        lightToIndex[scene.lights[i].get()] = i;

    // This point is not used in Lookup() for
    // PowerLightDistribution (which is the default strategy)
    const Point3f unusedPoint(0.f, 0.f, 0.f);
    const Distribution1D *lightDistr = lightDistribution->Lookup(unusedPoint);

    // Partition the image into tiles
    Film *film = camera->film;
    const Bounds2i sampleBounds = film->GetSampleBounds();
    const Vector2i sampleExtent = sampleBounds.Diagonal();
    const int tileSize = 16;
    const int nXTiles = (sampleExtent.x + tileSize - 1) / tileSize;
    const int nYTiles = (sampleExtent.y + tileSize - 1) / tileSize;

    // Record total number of pixels in the sample bound
    // This is the number of light subpaths and camera subpaths to be traced
    const int nPixels = sampleBounds.Area();

    // Prepare memory arena for each thread
    std::vector<MemoryArena> perThreadArenas(MaxThreadIndex());
    std::vector<MemoryArena> perThreadGridArenas(MaxThreadIndex());
    std::vector<MemoryArena> perThreadStage2Arenas(MaxThreadIndex());

    // Use number of samples per pixel as number of iterations
	int nIterations = sampler->samplesPerPixel;

    // Setup progress reporter
        ProgressReporter reporter((2 * nXTiles * nYTiles + 1) * nIterations,
                                  "Rendering");
	
    for (int iter = 0; iter < nIterations; iter++) {

        // ***** VCM STAGE 1 *****
        // Trace light subpaths and store light vertices to grid

        // ===== VCM Stage 1a: Trace light subpaths
		
        // Set up storage for light vertices
        std::vector<VCMVertex> *lightSubpaths = new std::vector<VCMVertex>[nPixels];

        if (scene.lights.size() > 0) {
            ParallelFor2D(
                [&](const Point2i tile) {
                    MemoryArena &arena = perThreadArenas[ThreadIndex];
                    int tileIndex = tile.y * nXTiles + tile.x;
                    std::unique_ptr<Sampler> tileSampler = sampler->Clone(tileIndex);

                    // Compute tile bounds
                    int x0 = sampleBounds.pMin.x + tile.x * tileSize;
                    int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
                    int y0 = sampleBounds.pMin.y + tile.y * tileSize;
                    int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
                    Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));

                    for (Point2i pPixel : tileBounds) {
                        tileSampler->StartPixel(pPixel);
                        tileSampler->SetSampleNumber(iter);

                        if (!InsideExclusive(pPixel, pixelBounds)) continue;
                        
                        VCMVertex *lightVertices =
                            arena.Alloc<VCMVertex>(maxDepth + 1);

                        // Linear interpolate light time
                        Float lightTime =
                            Lerp(tileSampler->Get1D(), camera->shutterOpen,
                                 camera->shutterClose);

                        // Trace light subpath for this pixel
                        int nLightVertices = VCMGenerateLightSubpath(
                            scene, *tileSampler, arena, vcmConstant,
                            maxDepth + 1, lightTime, *lightDistr, lightToIndex,
                            lightVertices);

                        // Put traced light vertices into storage vector
                        int pathIndex = pPixel.y * sampleExtent.x + pPixel.x;
                        for (int i = 0; i < nLightVertices; i++) {
                            lightVertices[i].lightSubPathIdx = pathIndex;
                            lightVertices[i].vertexInSubpathIndex = i;
                            lightSubpaths[pathIndex].push_back(
                                lightVertices[i]);
                        }
                    }

                    // Report completed a tile
                    reporter.Update();
                },
                Point2i(nXTiles, nYTiles));
        }

        // ===== End of Stage 1a
        
        // Note:
        // Stage 1b: connect light vertices to eye and
        // "Connect to a light source" in Stage 2a
        // do NOT need to be handled separately because ConnectBDPT() 
        // in the original bdpt implementation already takes these 
        // cases into account
        // We will handle these two cases in the same way

        // ===== VCM Stage 1c: Build hash grid for light vertices

        // Allocate voxel grid for storing light vertices
        const int hashSize = nPixels;
        std::vector<std::atomic<VCMVertexListNode *>> grid(hashSize);
        
        // Compute grid bounds
        int gridRes[3];
        Bounds3f gridBounds;
        for (int i = 0; i < nPixels; i++) {
            const std::vector<VCMVertex> lightSubpath = lightSubpaths[i];
            for (int s = 0; s < lightSubpath.size(); s++) {
                if (lightSubpath[s].beta.IsBlack()) continue;
                Bounds3f lightVertexBound =
                    Expand(Bounds3f(lightSubpath[s].p()), searchRadius);
                gridBounds = Union(gridBounds, lightVertexBound);
            }
        }

        // Compute resolution of grid in each dimension
        Vector3f diag = gridBounds.Diagonal();
        Float maxDiag = MaxComponent(diag);
        int baseGridRes = (int)(maxDiag / searchRadius);
        CHECK_GT(baseGridRes, 0);
        for (int i = 0; i < 3; i++) {
            gridRes[i] = std::max((int)(baseGridRes * diag[i] / maxDiag), 1);
        }

        // Add light vertices to grid
        ParallelFor(
            [&](int pathIndex) {
                MemoryArena &arena = perThreadGridArenas[ThreadIndex];
                std::vector<VCMVertex> &lightSubpath =
                    lightSubpaths[pathIndex];
                for (int s = 0; s < lightSubpath.size(); s++) {
                    if (!lightSubpath[s].beta.IsBlack()) {
                        VCMVertex &lightVertex = lightSubpath.at(s);
                        Point3f lightVertexPos = lightVertex.p();
                        Vector3f radiusExtend(searchRadius, searchRadius,
                                              searchRadius);

                        Point3i pMin, pMax;
                        ToGrid(lightVertexPos - radiusExtend, gridBounds,
                               gridRes, &pMin);
                        ToGrid(lightVertexPos + radiusExtend, gridBounds,
                               gridRes, &pMax);

                        for (int z = pMin.z; z <= pMax.z; ++z) {
                            for (int y = pMin.y; y <= pMax.y; ++y) {
                                for (int x = pMin.x; x <= pMax.x; ++x) {
                                    // Get grid hash
                                    int h = hash(Point3i(x, y, z), hashSize);

                                    // Create list node
                                    VCMVertexListNode *node =
                                        arena.Alloc<VCMVertexListNode>();
                                    node->lightVertex = &lightVertex;
                                    node->next = grid[h];

                                    // Atomically add node to the grid
                                    while (grid[h].compare_exchange_weak(
                                               node->next, node) == false)
                                        ;
                                }
                            }
                        }

                    }
                }
            },
            nPixels, 4096);

        // Report finished putting light vertices into grid
        reporter.Update();

        // ===== End of Stage 1c

        // ***** End of VCM STAGE 1 *****



        // ***** VCM STAGE 2 *****
        // Trace camera subpaths and handle VC and VM
        ParallelFor2D(
            [&](const Point2i tile) {
                MemoryArena &arena = perThreadStage2Arenas[ThreadIndex];
                int tileIndex = tile.y * nXTiles + tile.x;
                std::unique_ptr<Sampler> tileSampler =
                    sampler->Clone(2 * tileIndex);

                // Compute tile bounds
                int x0 = sampleBounds.pMin.x + tile.x * tileSize;
                int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
                int y0 = sampleBounds.pMin.y + tile.y * tileSize;
                int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
                Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));

                // Get film tile for the tile to add contribution later
                std::unique_ptr<FilmTile> filmTile =
                    camera->film->GetFilmTile(tileBounds);

                for (Point2i pPixel : tileBounds) {
                    tileSampler->StartPixel(pPixel);
                    tileSampler->SetSampleNumber(iter);

                    if (!InsideExclusive(pPixel, pixelBounds)) continue;

                    // Generate a film sample
                    Point2f pFilm = (Point2f)pPixel + tileSampler->Get2D();

                    VCMVertex *cameraVertices =
                        arena.Alloc<VCMVertex>(maxDepth + 2);

                    // Trace camera subpath for this pixel
                    int nCameraVertices = VCMGenerateCameraSubpath(
                        scene, *tileSampler, arena, vcmConstant, maxDepth + 1,
                        *camera, pFilm, cameraVertices);

                    // Get the corresponding light path for this pixel
                    const int pathIndex = pPixel.y * sampleExtent.x + pPixel.x;
                    std::vector<VCMVertex> lightSubpath = lightSubpaths[pathIndex];
                    VCMVertex *lightVertices = lightSubpath.data();
                    const int nLightVertices = lightSubpath.size();

                    Spectrum L(0.f);
                    for (int t = 1; t <= nCameraVertices; t++) {

                        // ===== VCM Stage 2a: Vertex Connection

                        for (int s = 0; s <= nLightVertices; s++) {
                            int depth = t + s - 2;
                            if ((s == 1 && t == 1) || depth < 0 ||
                                depth > maxDepth)
                                continue;

                            Point2f pFilmNew = pFilm;
                            Float vcMisWeight = 0.f;
                            VCMPathStrategy pathStrategy(VC);
                            Spectrum Lpath = VertexConnectionAndMerging(
                                scene, vcmConstant, pathStrategy, lightVertices,
                                cameraVertices, s, t, *lightDistr, lightToIndex,
                                *camera, *tileSampler, &pFilmNew, &vcMisWeight);

                            if (t != 1)
                                L += Lpath;
                            else
                                film->AddSplat(pFilmNew, Lpath);
                        }

                        // ===== End of Stage 2a

                        // ===== VCM Stage 2b: Vertex Merging

						if (t >= 2) {
                            // Check if camera vertex lies in grid
                            Point3i cameraVertexGridIndex;
                            if (ToGrid(cameraVertices[t - 1].p(), gridBounds,
                                       gridRes, &cameraVertexGridIndex)) {
                                // Get grid hash
                                int h = hash(cameraVertexGridIndex, hashSize);

                                // Get list of light vertices
                                for (VCMVertexListNode *node = grid[h].load(
                                         std::memory_order_relaxed);
                                     node != nullptr; node = node->next) {
                                    VCMVertex *lightVertex = node->lightVertex;

                                    // Skip if fewer than two vertices in light path
                                    if (lightVertex->vertexInSubpathIndex + 1 < 2)
                                        continue;

                                    // Check if light vertex is close enough
                                    if (DistanceSquared(cameraVertices[t - 1].p(),
                                                        lightVertex->p()) >
                                        searchRadius * searchRadius)
                                        continue;

									if (t + lightVertex->vertexInSubpathIndex + 1 >
                                        maxDepth)
                                        continue;

                                    // Get light path;
                                    VCMVertex *lightVertices =
                                        lightSubpaths[lightVertex->lightSubPathIdx].data();
                                    int s = lightVertex->vertexInSubpathIndex + 1;

                                    // Merge camera and light vertices
                                    Point2f unused = pFilm;
                                    Float vmMisWeight = 0.f;
                                    VCMPathStrategy pathStrategy(VM);
                                    Spectrum Lpath = VertexConnectionAndMerging(
                                        scene, vcmConstant, pathStrategy,
                                        lightVertices, cameraVertices, s, t,
                                        *lightDistr, lightToIndex, *camera,
                                        *tileSampler, &unused, &vmMisWeight);

                                    L += Lpath;
                                }
                            }
                        }
                        // ===== End of Stage 2b

                    } // Finished handling this pixel's paths
                    
                    filmTile->AddSample(pFilm, L);
                } // End of tile


                // Merge file tile into the main film
                film->MergeFilmTile(std::move(filmTile));

                // Report completed a tile
                reporter.Update();
            },
            Point2i(nXTiles, nYTiles));

        // ***** End of VCM STAGE 2 *****

        // Reset memory arenas for next iteration
        for (int i = 0; i < MaxThreadIndex(); i++) {
            perThreadArenas[i].Reset();
            perThreadGridArenas[i].Reset();
            perThreadStage2Arenas[i].Reset();
        }

        delete [] lightSubpaths;

    } // End of all iterations
    reporter.Done();
    film->WriteImage(1.0f / sampler->samplesPerPixel);
}

// [adopted from pbrt and modified]
Spectrum VertexConnectionAndMerging(
    const Scene &scene, VCMConstant vcmConstant, VCMPathStrategy pathStrategy,
    VCMVertex *lightVertices, VCMVertex *cameraVertices, int s, int t,
    const Distribution1D &lightDistr,
    const std::unordered_map<const Light *, size_t> &lightToIndex,
    const Camera &camera, Sampler &sampler, Point2f *pRaster,
    Float *misWeightPtr) {
    ProfilePhase _(Prof::BDPTConnectSubpaths);
    Spectrum L(0.f);
    // Ignore invalid connections related to infinite area lights
    if (t > 1 && s != 0 && cameraVertices[t - 1].type == VCMVertexType::Light)
        return Spectrum(0.f);

    // Perform connection and write contribution to _L_
    VCMVertex sampled;
    if (s == 0) {
        // Interpret the camera subpath as a complete path
        const VCMVertex &pt = cameraVertices[t - 1];
        if (pt.IsLight()) L = pt.Le(scene, cameraVertices[t - 2]) * pt.beta;
        DCHECK(!L.HasNaNs());
    } else if (t == 1) {
        // Sample a point on the camera and connect it to the light subpath
        const VCMVertex &qs = lightVertices[s - 1];
        if (qs.IsConnectible()) {
            VisibilityTester vis;
            Vector3f wi;
            Float pdf;
            Spectrum Wi = camera.Sample_Wi(qs.GetInteraction(), sampler.Get2D(),
                                           &wi, &pdf, pRaster, &vis);
            if (pdf > 0 && !Wi.IsBlack()) {
                // Initialize dynamically sampled vertex and _L_ for $t=1$ case
                sampled = VCMVertex::CreateCamera(&camera, vis.P1(), Wi / pdf);
                L = qs.beta * qs.f(sampled, TransportMode::Importance) * sampled.beta;
                if (qs.IsOnSurface()) L *= AbsDot(wi, qs.ns());
                DCHECK(!L.HasNaNs());
                // Only check visibility after we know that the path would
                // make a non-zero contribution.
                if (!L.IsBlack()) L *= vis.Tr(scene, sampler);
            }
        }
    } else if (s == 1) {
        // Sample a point on a light and connect it to the camera subpath
        const VCMVertex &pt = cameraVertices[t - 1];
        if (pt.IsConnectible()) {
            Float lightPdf;
            VisibilityTester vis;
            Vector3f wi;
            Float pdf;
            int lightNum =
                lightDistr.SampleDiscrete(sampler.Get1D(), &lightPdf);
            const std::shared_ptr<Light> &light = scene.lights[lightNum];
            Spectrum lightWeight = light->Sample_Li(
                pt.GetInteraction(), sampler.Get2D(), &wi, &pdf, &vis);
            if (pdf > 0 && !lightWeight.IsBlack()) {
                VCMEndpointInteraction ei(vis.P1(), light.get());
                sampled =
                    VCMVertex::CreateLight(ei, lightWeight / (pdf * lightPdf), 0);
                sampled.pdfFwd =
                    sampled.PdfLightOrigin(scene, pt, lightDistr, lightToIndex);
                L = pt.beta * pt.f(sampled, TransportMode::Radiance) * sampled.beta;
                if (pt.IsOnSurface()) L *= AbsDot(wi, pt.ns());
                // Only check visibility if the path would carry radiance.
                if (!L.IsBlack()) L *= vis.Tr(scene, sampler);
            }
        }
    } else {
        // Handle all other bidirectional connection cases
        const VCMVertex &qs = lightVertices[s - 1], &pt = cameraVertices[t - 1];
        if (qs.IsConnectible() && pt.IsConnectible()) {
            L = qs.beta * qs.f(pt, TransportMode::Importance) * pt.f(qs, TransportMode::Radiance) * pt.beta;
            //VLOG(2) << "General connect s: " << s << ", t: " << t <<
            //    " qs: " << qs << ", pt: " << pt << ", qs.f(pt): " << qs.f(pt, TransportMode::Importance) <<
            //    ", pt.f(qs): " << pt.f(qs, TransportMode::Radiance) << ", G: " << VCMG(scene, sampler, qs, pt) <<
            //    ", dist^2: " << DistanceSquared(qs.p(), pt.p());
            if (!L.IsBlack()) L *= VCMG(scene, sampler, qs, pt);
        }
    }

    ++totalPaths;
    if (L.IsBlack()) ++zeroRadiancePaths;
    ReportValue(pathLength, s + t - 2);

    // Compute MIS weight for connection strategy
    Float misWeight = L.IsBlack()
                          ? 0.f
                          : VCMMISWeight(scene, vcmConstant, pathStrategy,
                                         lightVertices, cameraVertices, sampled,
                                         s, t, lightDistr, lightToIndex);
    //VLOG(2) << "MIS weight for (s,t) = (" << s << ", " << t << ") connection: "
    //        << misWeight;
    DCHECK(!std::isnan(misWeight));
    L *= misWeight;
    if (misWeightPtr) *misWeightPtr = misWeight;
    return L;
}

// [adopted from pbrt and modified]
VCMIntegrator *CreateVCMIntegrator(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    Float searchRadius = params.FindOneFloat("searchRadius", 1.f);
    bool visualizeStrategies = params.FindOneBool("visualizestrategies", false);
    bool visualizeWeights = params.FindOneBool("visualizeweights", false);

    if ((visualizeStrategies || visualizeWeights) && maxDepth > 5) {
        Warning(
            "visualizestrategies/visualizeweights was enabled, limiting "
            "maxdepth to 5");
        maxDepth = 5;
    }
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }

    std::string lightStrategy = params.FindOneString("lightsamplestrategy",
                                                     "power");
    return new VCMIntegrator(sampler, camera, maxDepth, searchRadius,
                             visualizeStrategies, visualizeWeights, pixelBounds,
                             lightStrategy);
}

}  // namespace pbrt
