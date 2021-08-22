PatchWarp relies on the following publicly available codes.

## Akinori_Mitani_motion_correct_files
Pyramid registration approach for rigid motion corrections by Akinori Mitani. The full original codes are availabe at https://github.com/amitani/matlab_motion_correct. The details can be found in the following paper: Mitani & Komiyama, "Real-Time Processing of Two-Photon Calcium Imaging Data Including Lateral Motion Artifact Correction", Front. Neuroinform., 18, 2018 

## ecc_patchwarp
ECC is a gradient-based image registration algorithm. PatchWarp uses ECC to derive affine transformation matrices. The original codes are available at https://www.mathworks.com/matlabcentral/fileexchange/27253. Some functions were borrowed from its variant https://www.mathworks.com/matlabcentral/fileexchange/62921. The details can be found in the following paper: G.D.Evangelidis, E.Z.Psarakis, "Parametric Image Alignment using Enhanced Correlation Coefficient", IEEE Trans. on PAMI, vol.30, no.10, 2008"

## ScanImageTiffReader
Fast tiff readers.
https://vidriotech.gitlab.io/scanimagetiffreader-matlab/
