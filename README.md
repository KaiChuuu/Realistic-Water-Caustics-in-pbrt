# Realistic Water Caustics in pbrt
Completed a grad research project with teammate Ivan.
 
The goal of the project was to render physically accurate caustic effects through bodies of water in pbrt.

We added our work ontop of <a href="https://pbrt.org/">pbrt-v3</a>, utilizing their Vertex Connection in BDPT and Vertex Merging in PM to recreate <a href="http://iliyan.com/publications/ImplementingVCM">Vertex Connection and Merging (VCM)</a>.

## What I contributed

I worked on the stage 1 (a and c steps) for VCM (In <a href="/VCM integration into pbrt/vcm.cpp">vcm.cpp</a> from line 530).  
This involved tracing light subpaths from all light sources, creating a vector grid and storing the light vertices within the grid to be used for stage 2.

<img src="/Readme Images/VCM overview.jpg" alt="VCM Overview" width="60%">
<i>Obtained from <a href="http://iliyan.com/publications/ImplementingVCM">VCM</a> paper</i>

<br/><br/>
I also worked on creating all the scenes in blender and,  
tested our vcm implementation against other rendering methods in various lighting and scene conditions.

<img src="/Readme Images/Rendered Test Image Compilation.png" alt="Rendered Test Image Compilation" width="80%">
<i>Compilation of some of the rendered scenes during testing.</i>

<br/><br/>
## File Navigation

Git contains file presentation, report and submission code for implmentation of VCM into pbrt.

In more detail... (from top to bottom)
 * Final rendered images using our implemented vcm method.
 * Contains pbrt generated files from our own created scenes and vcm method.
 * Test scenes created in blender, plus some pbrt generated test scenes.  
   (Blender scenes do not contain water fluid data, so it will need to be reconstructed if repurposed.)
 * vcm and api file integrated into pbrt.

For a general understanding of the project,  
The <a href="https://github.com/RedDogClifford/Realistic-Water-Caustics-in-pbrt/blob/main/CMPT985_Project_Final_Presentation.pdf">final project presentation</a> contains a very good summary of all the work done for the grad project.

