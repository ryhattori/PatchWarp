![logo](https://user-images.githubusercontent.com/25396523/130375407-d5a7646c-3b4e-42cb-baa5-268f02f68595.png)


PatchWarp is an image processing pipeline for neuronal calcium imaging data. It can correct complex image distortions that slowly occur during long imaging sessions. First, the pipeline performs rigid motion corrections by iterative re-estimation of template images. Then, the imaging field is split into user-specified number of subfields. A gradient-based algorithm independently finds the best affine transformation matrix for each frame of each subfield to correct the across-time image distortion of each subfield. The distortion-corrected subfields are stitched together like patchwork to reconstruct the distortion-corrected whole imaging field. PatchWarp can be also used to register images from different imaging sessions for longitudinal activity analyses.

---
**Before (Left) and after (Right) PatchWarp warp correction for within-session image distortions**    
<img src="https://user-images.githubusercontent.com/25396523/131230196-1938d133-6ea5-4814-af53-41e9a949ddae.gif" width="720" height="360">  
(2.25 hrs *in vivo* 2-photon calcium imaging of cell bodies with complex distortions)  
<img src="https://user-images.githubusercontent.com/25396523/155606780-a26a148f-3cb9-4083-97b1-3d28fe36f465.gif" width="720" height="360">  
(19 min *in vivo* 2-photon calcium imaging of axons with complex distortions. Image credit to Chi Ren.)  

**Before (Left) and after (Right) PatchWarp across-session image registration**   
<img src="https://user-images.githubusercontent.com/25396523/134836357-30dc6772-b6a7-487e-83b5-adc272076db9.jpg" width="360" height="360"> <img src="https://user-images.githubusercontent.com/25396523/134836358-9c686950-db99-45e5-a43f-d596964c09bf.jpg" width="360" height="360">  
(A later imaging session (cyan) was registered to an earlier imaging session (red)) 

## Installation
Download files from this github repository, and add all files to your MATLAB path. 

## How to use
Please check example demo files. 

For within-session distortion correction, please check [**patchwarp_demo.m**](https://github.com/ryhattori/PatchWarp/blob/main/PatchWarp/patchwarp_demo.m).   
Example data to run the demo code is available at https://doi.org/10.6084/m9.figshare.19212501.v1

For across-session image registration, please check [**patchwarp_across_sessions_demo.m**](https://github.com/ryhattori/PatchWarp/blob/main/PatchWarp/patchwarp_across_sessions_demo.m).  
Example data to run the demo code is available at https://doi.org/10.6084/m9.figshare.19217421.v1

## Citation
Example citation format for the preprint:  
- Hattori, R. and Komiyama, T. PatchWarp: Corrections of non-uniform image distortions in two-photon calcium imaging data by patchwork affine transformations. **_Cell Reports Methods_** (2022), https://doi.org/10.1016/j.crmeth.2022.100205.

Software DOI from Zenodo:   
- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5232757.svg)](https://doi.org/10.5281/zenodo.5232757)
