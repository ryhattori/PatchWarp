# PatchWarp
PatchWarp is an image processing pipeline for neuronal calcium imaging data. It can correct complex image distortions that slowly occur during long imaging sessions. First, the pipeline performs rigid motion corrections by iterative re-estimation of template images. Then, the imaging field is split into user-specified number of subfields. A gradient descent algorithm independently finds the best affine transformation matrix for each frame of each subfield to correct the across-time image distortion of each subfield. The distortion-corrected subfields are stitched together like patchwork to reconstruct the distortion-corrected whole imaging field.


https://user-images.githubusercontent.com/25396523/130366052-f2bb202c-9b83-4884-ae0c-0203ffee81a3.mp4


https://user-images.githubusercontent.com/25396523/130366460-96bd5d68-a089-489b-8bce-eb7948d63a34.mov




## Installation
Download and add all files to MATLAB path.

## How to use
Please check an example demo in **patchwarp_demo.mat** file.

## Citation
Example citation format:  
Hattori, R. (2021). PatchWarp (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.5232758

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5232758.svg)](https://doi.org/10.5281/zenodo.5232758)
