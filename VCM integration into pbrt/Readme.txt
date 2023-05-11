Please follow the steps below to compile and run pbrt with our VCM module.



1. Clone the official pbrt repo from https://github.com/mmp/pbrt-v3. Be sure to use the --recursive flag to download all dependencies as well.

2. Put our implementation files in the following locations:
   a) vcm.cpp and vcm.h at \pbrt-v3\src\integrators
   b) api.cpp at \pbrt-v3\src\core (replace the existing file)

3. Follow the instructions at https://github.com/mmp/pbrt-v3 to generate a solution project for your machine using CMake.

4. Continue following the instructions at https://github.com/mmp/pbrt-v3 to compile the project.

5. You can now run pbrt with our VCM module. You can render an image with:
   .\pbrt.exe <path to a scene's .pbrt file>

6. To switch a scene to use our VCM integrator, open the scene's .pbrt file and replace the existing integrator with:
   Integrator "vcm"
   
7. You can specific the search radius for vertex merging and max depth for sub path generation like this:
   Integrator "vcm" "float searchRadius" [1.0] "integer maxdepth" [5]

8. You can also specific the number of iterations vcm runs by setting the number of pixel samples at the sampler. For example:
   Sampler "sobol" "integer pixelsamples" [16]

For api.cpp, the only modification is in MakeIntegrator().
For vcm.cpp/.h, please refer to the comments in the code (especially before each method's declaration) to see the amount of work we did on each method.

*This submission package does not include the blender files for the scenes we created as they include fluid simulation files and are over 5 gigs in size.
They are available upon request.

Please don't hesitate to reach out to us if you have any questions.