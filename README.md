![logo](https://user-images.githubusercontent.com/25396523/130375407-d5a7646c-3b4e-42cb-baa5-268f02f68595.png)


PatchWarp is an image processing pipeline for neuronal calcium imaging data. It can correct complex image distortions that slowly occur during long imaging sessions. First, the pipeline performs rigid motion corrections by iterative re-estimation of template images. Then, the imaging field is split into user-specified number of subfields. A gradient-based algorithm independently finds the best affine transformation matrix for each frame of each subfield to correct the across-time image distortion of each subfield. The distortion-corrected subfields are stitched together like patchwork to reconstruct the distortion-corrected whole imaging field.

---
**Before PatchWarp affine transformations (Left)**    /    **After PatchWarp affine transformations (Right)**
<img src="https://user-images.githubusercontent.com/25396523/130368392-0e2c448c-7a9c-4458-9a73-20d63ca06694.gif" width="680" height="340">  
(Example 2.25 hrs *in vivo* two-photon calcium imaging with complex image distortions)

## Installation
Download and add all files to MATLAB path.

## How to use
Please check an example demo in [**patchwarp_demo.mat**](https://github.com/ryhattori/PatchWarp/blob/main/PatchWarp/patchwarp_demo.m) file.

## Citation
Example citation format:  
- Hattori, R. (2021). PatchWarp (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.5232758

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5232758.svg)](https://doi.org/10.5281/zenodo.5232758)

This pipeline was first used for the image processings for the following paper:  
- Hattori, R., Danskin, B., Babic, Z., Mlynaryk, N., and Komiyama, T.  
Area-Specificity and Plasticity of History-Dependent Value Coding During Learning.  
Cell, 2019 Jun 13;177(7):1858-1872.e15. doi: 10.1016/j.cell.2019.04.027.  
[Link to the publication](https://www.cell.com/cell/fulltext/S0092-8674(19)30446-5)
