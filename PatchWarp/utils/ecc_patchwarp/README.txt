This software implements the ECC image alignment algorithm as it is
presented in the paper 

G.D.Evangelidis, E.Z.Psarakis, "Parametric Image Alignment
using Enhanced Correlation Coefficient", IEEE Trans. on PAMI, vol.30, no.10, 2008"

Given two images (template, input), ECC algorithm computes the 
transformation (translation, euclidean, affine or homography) between the images, 
so that, when it is applied to the input 
image, provides a new warped image close to the template one. 

The code is provided for research proposes only. If you use this software, 
please cite the above paper.

--------------------------------------------------------------------------------
Usage:
To see how to use the code, please see the file ecc_demo or type 'help ecc'.
For a detailed review, see the above paper.

Note that the algorithm ignores a border of the template whose size is equal to
five percent of M, where M = mean([size(template,1), size(template,2)]).
You might want to include the whole template area, so edit the ecc.m
file and set the "margin" variable equal to zero.
