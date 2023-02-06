# optoacoustic_SAFT

This code implements a simulation environment for optoacoustic microspcopy scans for concave focused transducers. Scans can be performed on arbitrary objects
by providing an image representing the scanned object.\
Synthetic Aperture Fcousing Technique (SAFT) algorithms are implemented as well. As the purpose of this work was to test the limitations of SAFT algorithms. 

## Important Files
**Presentations/SAFT Final Report** - This is a document explaining the methods used in the code as well as important concepts regarding SAFT. The document also shows the main results of our work in 2D and 3D showing SAFT quality degradation when imaging objects angled in the depth direction.

## Run Examples
**Code/Test_3D.m** - a code example for scanning depth angled lines and applying SAFT in the x direction, y direction and both in a cross pattern. The test is configured for imaging a long angled line and the angle and length can be adjusted to recreate all the results.\
**Code/Test_2D.m** - a code example for scanning depth angled lines and applying SAFT in the y direction using different SAFT methods: Basic, SIR, spherical SIR, and Coherence factor (CF).

## Scans 
**Sinogram.m** - implements a 2D scan of a given image.\
**scan_3D.m** - implements a 3D scan of a given image, by translating the image and the transducer's response function.

## SAFT methods 
All methods are inside the folder Code/SAFT_Methods.\
Methods that contain 3D in the name are for 3D scans and the rest are for 2D.\
Each method is briefly explained in the documentation and broadly explained in the report document.
