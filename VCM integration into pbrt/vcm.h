
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_VCM_H
#define PBRT_INTEGRATORS_VCM_H

// integrators/vcm.h*
#include <unordered_map>
#include "camera.h"
#include "integrator.h"
#include "interaction.h"
#include "light.h"
#include "pbrt.h"
#include "reflection.h"
#include "sampling.h"
#include "scene.h"

namespace pbrt {

/// Forward declaration (correction term for adjoint BSDF with shading normals)
extern Float VCMCorrectShadingNormal(const SurfaceInteraction &isect,
                                  const Vector3f &wo, const Vector3f &wi,
                                  TransportMode mode);

// VCMEndpointInteraction Declarations
// [Directly adopted from pbrt]
struct VCMEndpointInteraction : Interaction {
    union {
        const Camera *camera;
        const Light *light;
    };
    // VCMEndpointInteraction Public Methods
    VCMEndpointInteraction() : Interaction(), light(nullptr) {}
    VCMEndpointInteraction(const Interaction &it, const Camera *camera)
        : Interaction(it), camera(camera) {}
    VCMEndpointInteraction(const Camera *camera, const Ray &ray)
        : Interaction(ray.o, ray.time, ray.medium), camera(camera) {}
    VCMEndpointInteraction(const Light *light, const Ray &r, const Normal3f &nl)
        : Interaction(r.o, r.time, r.medium), light(light) {
        n = nl;
    }
    VCMEndpointInteraction(const Interaction &it, const Light *light)
        : Interaction(it), light(light) {}
    VCMEndpointInteraction(const Ray &ray)
        : Interaction(ray(1), ray.time, ray.medium), light(nullptr) {
        n = Normal3f(-ray.d);
    }
};

// VCM Helper Definitions
enum class VCMVertexType { Camera, Light, Surface, Medium };
struct VCMVertex;

// [Directly adopted from pbrt]
template <typename Type>
class VCMScopedAssignment {
  public:
    // VCMScopedAssignment Public Methods
    VCMScopedAssignment(Type *target = nullptr, Type value = Type())
        : target(target) {
        if (target) {
            backup = *target;
            *target = value;
        }
    }
    ~VCMScopedAssignment() {
        if (target) *target = backup;
    }
    VCMScopedAssignment(const VCMScopedAssignment &) = delete;
    VCMScopedAssignment &operator=(const VCMScopedAssignment &) = delete;
    VCMScopedAssignment &operator=(VCMScopedAssignment &&other) {
        if (target) *target = backup;
        target = other.target;
        backup = other.backup;
        other.target = nullptr;
        return *this;
    }

  private:
    Type *target, backup;
};

// [Directly adopted from pbrt]
inline Float VCMInfiniteLightDensity(
    const Scene &scene, const Distribution1D &lightDistr,
    const std::unordered_map<const Light *, size_t> &lightToDistrIndex,
    const Vector3f &w) {
    Float pdf = 0;
    for (const auto &light : scene.infiniteLights) {
        CHECK(lightToDistrIndex.find(light.get()) != lightToDistrIndex.end());
        size_t index = lightToDistrIndex.find(light.get())->second;
        pdf += light->Pdf_Li(Interaction(), -w) * lightDistr.func[index];
    }
    return pdf / (lightDistr.funcInt * lightDistr.Count());
}

// [Completely New]
struct VCMConstant {
    Float nVC;
    Float nVM;
    Float nVCNormalizationFactor;
    Float nVMNormalizationFactor;
    Float etaVCM;
    Float invEtaVCM;
    int nLight;
};

// VCM Declarations
// [Adopted from pbrt and substantially modified]
class VCMIntegrator : public Integrator {
  public:
    // VCMIntegrator Public Methods
    VCMIntegrator(std::shared_ptr<Sampler> sampler,
                   std::shared_ptr<const Camera> camera,
                   int maxDepth, Float searchRadius,
                   bool visualizeStrategies, bool visualizeWeights,
                   const Bounds2i &pixelBounds,
                   const std::string &lightSampleStrategy = "power")
        : sampler(sampler),
          camera(camera),
          maxDepth(maxDepth),
          searchRadius(searchRadius),
          visualizeStrategies(visualizeStrategies),
          visualizeWeights(visualizeWeights),
          pixelBounds(pixelBounds),
          lightSampleStrategy(lightSampleStrategy) {

        vcmConstant.nVC = 1;
        vcmConstant.nVM = camera->film->GetSampleBounds().Area();
        vcmConstant.nVCNormalizationFactor = 1.f / vcmConstant.nVC;
        vcmConstant.nVMNormalizationFactor = 1.f / vcmConstant.nVM;
        vcmConstant.etaVCM = vcmConstant.nVM / vcmConstant.nVC * Pi *
                              searchRadius * searchRadius;
        vcmConstant.invEtaVCM = 1.f / vcmConstant.etaVCM;
        vcmConstant.nLight = camera->film->GetSampleBounds().Area();
    }
    void Render(const Scene &scene);

  private:
    // VCMIntegrator Private Data
    std::shared_ptr<Sampler> sampler;
    std::shared_ptr<const Camera> camera;
    const int maxDepth;
    const Float searchRadius;
    const bool visualizeStrategies;
    const bool visualizeWeights;
    const Bounds2i pixelBounds;
    const std::string lightSampleStrategy;
    VCMConstant vcmConstant;
};

// [Completely New]
struct VCMVertexListNode {
    VCMVertex *lightVertex;
    VCMVertexListNode *next;
};

// [Adopted from pbrt and modified]
struct VCMVertex {
    // VCMVertex Public Data
    VCMVertexType type;
    Spectrum beta;
#ifdef PBRT_HAVE_NONPOD_IN_UNIONS
    union {
#else
    struct {
#endif  // PBRT_HAVE_NONPOD_IN_UNIONS
        VCMEndpointInteraction ei;
        MediumInteraction mi;
        SurfaceInteraction si;
    };
    bool delta = false;
    Float pdfPos = 0, pdfDir = 0;

	// [Completely New]
    int lightSubPathIdx = 0;
    int vertexInSubpathIndex = 0;
    Float pdfFwd = 0, pdfRev = 0;
    Float pdfSelectLight = 0;
    Float dVCM = 0, dVC = 0, dVM = 0;

    // VCMVertex Public Methods
    VCMVertex() : ei() {}
    VCMVertex(VCMVertexType type, const VCMEndpointInteraction &ei, const Spectrum &beta)
        : type(type), beta(beta), ei(ei) {}
    VCMVertex(const SurfaceInteraction &si, const Spectrum &beta)
        : type(VCMVertexType::Surface), beta(beta), si(si) {}

    // Need to define these two to make compilers happy with the non-POD
    // objects in the anonymous union above.
    VCMVertex(const VCMVertex &v) { memcpy(this, &v, sizeof(VCMVertex)); }
    VCMVertex &operator=(const VCMVertex &v) {
        memcpy(this, &v, sizeof(VCMVertex));
        return *this;
    }

    // Starting from a camera
    static inline VCMVertex CreateCamera(const Camera *camera, const Ray &ray,
                                         const Spectrum &beta,
                                         const Float pdfPos,
                                         const Float pdfDir);
    // Ending in a camera
    static inline VCMVertex CreateCamera(const Camera *camera,
                                         const Interaction &it,
                                         const Spectrum &beta);
    // Starting from a light source
    static inline VCMVertex CreateLight(const Light *light, const Ray &ray,
                                        const Normal3f &nLight,
                                        const Spectrum &Le, Float pdf,
                                        const Float pdfSelectLight,
                                        const Float pdfPos, const Float pdfDir);
    // Ending in a light source
    // Note: Infinite area light
    static inline VCMVertex CreateLight(const VCMEndpointInteraction &ei,
                                        const Spectrum &beta, Float pdf);
    static inline VCMVertex CreateMedium(const MediumInteraction &mi,
                                         const Spectrum &beta, Float pdf,
                                         const VCMVertex &prev);
    static inline VCMVertex CreateSurface(const SurfaceInteraction &si,
                                          const Spectrum &beta, Float pdf,
                                          const VCMVertex &prev);
    VCMVertex(const MediumInteraction &mi, const Spectrum &beta)
        : type(VCMVertexType::Medium), beta(beta), mi(mi) {}
    const Interaction &GetInteraction() const {
        switch (type) {
        case VCMVertexType::Medium:
            return mi;
        case VCMVertexType::Surface:
            return si;
        default:
            return ei;
        }
    }
    const Point3f &p() const { return GetInteraction().p; }
    Float time() const { return GetInteraction().time; }
    const Normal3f &ng() const { return GetInteraction().n; }
    const Normal3f &ns() const {
        if (type == VCMVertexType::Surface)
            return si.shading.n;
        else
            return GetInteraction().n;
    }
    bool IsOnSurface() const { return ng() != Normal3f(); }
    Spectrum f(const VCMVertex &next, TransportMode mode) const {
        Vector3f wi = next.p() - p();
        if (wi.LengthSquared() == 0) return 0.;
        wi = Normalize(wi);
        switch (type) {
        case VCMVertexType::Surface:
            return si.bsdf->f(si.wo, wi) *
                VCMCorrectShadingNormal(si, si.wo, wi, mode);
        case VCMVertexType::Medium:
            return mi.phase->p(mi.wo, wi);
        default:
            LOG(FATAL) << "VCMVertex::f(): Unimplemented";
            return Spectrum(0.f);
        }
    }
    bool IsConnectible() const {
        switch (type) {
        case VCMVertexType::Medium:
            return true;
        case VCMVertexType::Light:
            return (ei.light->flags & (int)LightFlags::DeltaDirection) == 0;
        case VCMVertexType::Camera:
            return true;
        case VCMVertexType::Surface:
            return si.bsdf->NumComponents(BxDFType(BSDF_DIFFUSE | BSDF_GLOSSY |
                                                   BSDF_REFLECTION |
                                                   BSDF_TRANSMISSION)) > 0;
        }
        LOG(FATAL) << "Unhandled vertex type in IsConnectable()";
        return false;  // NOTREACHED
    }
    bool IsLight() const {
        return type == VCMVertexType::Light ||
               (type == VCMVertexType::Surface && si.primitive->GetAreaLight());
    }
    bool IsDeltaLight() const {
        return type == VCMVertexType::Light && ei.light &&
               pbrt::IsDeltaLight(ei.light->flags);
    }
    bool IsInfiniteLight() const {
        return type == VCMVertexType::Light &&
               (!ei.light || ei.light->flags & (int)LightFlags::Infinite ||
                ei.light->flags & (int)LightFlags::DeltaDirection);
    }
    Spectrum Le(const Scene &scene, const VCMVertex &v) const {
        if (!IsLight()) return Spectrum(0.f);
        Vector3f w = v.p() - p();
        if (w.LengthSquared() == 0) return 0.;
        w = Normalize(w);
        if (IsInfiniteLight()) {
            // Return emitted radiance for infinite light sources
            Spectrum Le(0.f);
            for (const auto &light : scene.infiniteLights)
                Le += light->Le(Ray(p(), -w));
            return Le;
        } else {
            const AreaLight *light = si.primitive->GetAreaLight();
            CHECK(light != nullptr);
            return light->L(si, w);
        }
    }
    friend std::ostream &operator<<(std::ostream &os, const VCMVertex &v) {
        return os << v.ToString();
    }
    std::string ToString() const {
        std::string s = std::string("[VCMVertex type: ");
        switch (type) {
        case VCMVertexType::Camera:
            s += "camera";
            break;
        case VCMVertexType::Light:
            s += "light";
            break;
        case VCMVertexType::Surface:
            s += "surface";
            break;
        case VCMVertexType::Medium:
            s += "medium";
            break;
        }
        s += std::string(" connectible: ") +
            std::string(IsConnectible() ? "true" : "false");
        s += StringPrintf("\n  p: [ %f, %f, %f ] ng: [ %f, %f, %f ]", p().x, p().y,
                          p().z, ng().x, ng().y, ng().z);
        s += StringPrintf("\n  pdfFwd: %f pdfRev: %f beta: ", pdfFwd, pdfRev) +
             beta.ToString();
        switch (type) {
        case VCMVertexType::Camera:
            // TODO
            break;
        case VCMVertexType::Light:
            // TODO
            break;
        case VCMVertexType::Surface:
            s += std::string("\n  bsdf: ") + si.bsdf->ToString();
            break;
        case VCMVertexType::Medium:
            s += std::string("\n  phase: ") + mi.phase->ToString();
            break;
        }
        s += std::string(" ]");
        return s;
    }
    Float ConvertDensity(Float pdf, const VCMVertex &next) const {
        // Return solid angle density if _next_ is an infinite area light
        if (next.IsInfiniteLight()) return pdf;
        Vector3f w = next.p() - p();
        if (w.LengthSquared() == 0) return 0;
        Float invDist2 = 1 / w.LengthSquared();
        if (next.IsOnSurface())
            pdf *= AbsDot(next.ng(), w * std::sqrt(invDist2));
        return pdf * invDist2;
    }
    Float Pdf(const Scene &scene, const VCMVertex *prev,
              const VCMVertex &next, const bool area = true) const {
        if (type == VCMVertexType::Light) return PdfLight(scene, next);
        // Compute directions to preceding and next vertex
        Vector3f wn = next.p() - p();
        if (wn.LengthSquared() == 0) return 0;
        wn = Normalize(wn);
        Vector3f wp;
        if (prev) {
            wp = prev->p() - p();
            if (wp.LengthSquared() == 0) return 0;
            wp = Normalize(wp);
        } else
            CHECK(type == VCMVertexType::Camera);

        // Compute directional density depending on the vertex types
        Float pdf = 0, unused;
        if (type == VCMVertexType::Camera)
            ei.camera->Pdf_We(ei.SpawnRay(wn), &unused, &pdf);
        else if (type == VCMVertexType::Surface)
            pdf = si.bsdf->Pdf(wp, wn);
        else if (type == VCMVertexType::Medium)
            pdf = mi.phase->p(wp, wn);
        else
            LOG(FATAL) << "VCMVertex::Pdf(): Unimplemented";

        if (area) {
            // Return probability per unit area at vertex _next_
            return ConvertDensity(pdf, next);
        } else {
            return pdf;
        }
    }
    Float PdfLight(const Scene &scene, const VCMVertex &v) const {
        Vector3f w = v.p() - p();
        Float invDist2 = 1 / w.LengthSquared();
        w *= std::sqrt(invDist2);
        Float pdf;
        if (IsInfiniteLight()) {
            // Compute planar sampling density for infinite light sources
            Point3f worldCenter;
            Float worldRadius;
            scene.WorldBound().BoundingSphere(&worldCenter, &worldRadius);
            pdf = 1 / (Pi * worldRadius * worldRadius);
        } else {
            // Get pointer _light_ to the light source at the vertex
            CHECK(IsLight());
            const Light *light = type == VCMVertexType::Light
                                     ? ei.light
                                     : si.primitive->GetAreaLight();
            CHECK(light != nullptr);

            // Compute sampling density for non-infinite light sources
            Float pdfPos, pdfDir;
            light->Pdf_Le(Ray(p(), w, Infinity, time()), ng(), &pdfPos, &pdfDir);
            pdf = pdfDir * invDist2;
        }
        if (v.IsOnSurface()) pdf *= AbsDot(v.ng(), w);
        return pdf;
    }
    Float PdfLightOrigin(const Scene &scene, const VCMVertex &v,
                         const Distribution1D &lightDistr,
                         const std::unordered_map<const Light *, size_t>
                             &lightToDistrIndex) const {
        Vector3f w = v.p() - p();
        if (w.LengthSquared() == 0) return 0.;
        w = Normalize(w);
        if (IsInfiniteLight()) {
            // Return solid angle density for infinite light sources
            return VCMInfiniteLightDensity(scene, lightDistr, lightToDistrIndex,
                                        w);
        } else {
            // Return solid angle density for non-infinite light sources
            Float pdfPos, pdfDir, pdfChoice = 0;

            // Get pointer _light_ to the light source at the vertex
            CHECK(IsLight());
            const Light *light = type == VCMVertexType::Light
                                     ? ei.light
                                     : si.primitive->GetAreaLight();
            CHECK(light != nullptr);

            // Compute the discrete probability of sampling _light_, _pdfChoice_
            CHECK(lightToDistrIndex.find(light) != lightToDistrIndex.end());
            size_t index = lightToDistrIndex.find(light)->second;
            pdfChoice = lightDistr.DiscretePDF(index);

            light->Pdf_Le(Ray(p(), w, Infinity, time()), ng(), &pdfPos, &pdfDir);
            return pdfPos * pdfChoice;
        }
    }
};

// [Adopted from pbrt and modified]
extern int VCMGenerateCameraSubpath(const Scene &scene, Sampler &sampler,
                                    MemoryArena &arena, VCMConstant vcmConstant,
                                    int maxDepth, const Camera &camera,
                                    const Point2f &pFilm, VCMVertex *path);

// [Adopted from pbrt and modified]
extern int VCMGenerateLightSubpath(
    const Scene &scene, Sampler &sampler, MemoryArena &arena,
    VCMConstant vcmConstant, int maxDepth, Float time,
    const Distribution1D &lightDistr,
    const std::unordered_map<const Light *, size_t> &lightToIndex,
    VCMVertex *path);

// [Completely New]
static enum VCMPathStrategy {VC, VM};

// [Adopted from pbrt and modified]
Spectrum VertexConnectionAndMerging(
    const Scene &scene, VCMConstant vcmConstant, VCMPathStrategy pathStrategy,
    VCMVertex *lightVertices, VCMVertex *cameraVertices, int s, int t,
    const Distribution1D &lightDistr,
    const std::unordered_map<const Light *, size_t> &lightToIndex,
    const Camera &camera, Sampler &sampler, Point2f *pRaster,
    Float *misWeight = nullptr);

// [Adopted from pbrt and substantially modified]
Spectrum VCMVertexMerging(const Scene &scene, VCMVertex *cameraVertex,
                          VCMVertex *cameraVertexPrev, VCMVertex *lightVertex,
                          VCMVertex *lightVertexPrev, VCMConstant vcmConstants);

// [Adopted from pbrt and modified]
VCMIntegrator *CreateVCMIntegrator(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera);


// Everything below: [Adopted from pbrt and modified]
// VCMVertex Inline Method Definitions
// Starting from a camera
inline VCMVertex VCMVertex::CreateCamera(const Camera *camera, const Ray &ray,
                                         const Spectrum &beta,
                                         const Float pdfPos,
                                         const Float pdfDir) {
    VCMVertex v = VCMVertex(VCMVertexType::Camera, VCMEndpointInteraction(camera, ray), beta);
    v.pdfPos = pdfPos;
    v.pdfDir = pdfDir;
    return v;
}

// Ending in a camera
inline VCMVertex VCMVertex::CreateCamera(const Camera *camera, const Interaction &it,
                                   const Spectrum &beta) {
    VCMVertex v = VCMVertex(VCMVertexType::Camera, VCMEndpointInteraction(it, camera), beta);
    return v;
}

// Starting from a light source
inline VCMVertex VCMVertex::CreateLight(const Light *light, const Ray &ray,
                                        const Normal3f &Nl, const Spectrum &Le,
                                        Float pdf, const Float pdfSelectLight,
                                        const Float pdfPos,
                                        const Float pdfDir) {
    VCMVertex v(VCMVertexType::Light, VCMEndpointInteraction(light, ray, Nl), Le);
    v.pdfFwd = pdf;
    v.pdfSelectLight = pdfSelectLight;
    v.pdfPos = pdfPos;
    v.pdfDir = pdfDir;
    return v;
}

inline VCMVertex VCMVertex::CreateSurface(const SurfaceInteraction &si,
                                    const Spectrum &beta, Float pdf,
                                          const VCMVertex &prev) {
    VCMVertex v(si, beta);
    v.pdfFwd = prev.ConvertDensity(pdf, v);
    return v;
}

inline VCMVertex VCMVertex::CreateMedium(const MediumInteraction &mi,
                                   const Spectrum &beta, Float pdf,
                                         const VCMVertex &prev) {
    VCMVertex v(mi, beta);
    v.pdfFwd = prev.ConvertDensity(pdf, v);
    return v;
}

// Ending in a light source
// Note: Infinite area light
inline VCMVertex VCMVertex::CreateLight(const VCMEndpointInteraction &ei,
                                        const Spectrum &beta, Float pdf) {
    VCMVertex v(VCMVertexType::Light, ei, beta);
    v.pdfFwd = pdf;
    return v;
}

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_VCM_H
