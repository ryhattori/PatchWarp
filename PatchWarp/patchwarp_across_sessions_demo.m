%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PatchWarp application to registration between different imaging sessions
% -------------------
% Note that this across-session registration function is standalone. You do
% not need to process your imaging data using the main PatchWarp pipeline
% to use this across-session registration function. You only need 2D images
% such as the mean image or the max-projection image for each imaging session. 
% 
% Released by Ryoma Hattori
% Email: rhattori0204@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add PatchWarp directory to MATLAB path
patchwarp_path = 'Z:\People\Ryoma\PatchWarp';
addpath(genpath(patchwarp_path))

%% Specifiy paths of the summary images from 2 different sessions
% Summary images can be either the mean image of all frames (e.g. downsampled_50_warped_mean.tif) or 
% the max-projection of downsampled movie (e.g. downsampled_50_warped_mean.tif) from the outputs of PatchWarp pipeline.
% They do not need to be from the main PatchWarp pipeline if you used a different motion correction software.
% The results are more robust and accurate if both mean and max-projection images are used.
session1_mean_path = 'G:\Data\171024\RH825\post_warp_affine\downsampled\downsampled_50_warped_mean.tif';
session1_max_path = 'G:\Data\171024\RH825\post_warp_affine\downsampled\downsampled_50_warped_max.tif';
session2_mean_path = 'G:\Data\171026\RH825\post_warp_affine\downsampled\downsampled_50_warped_mean.tif';
session2_max_path = 'G:\Data\171026\RH825\post_warp_affine\downsampled\downsampled_50_warped_max.tif';

%% Specify saving directory
save_path = 'G:\Data\171026\RH825';    % Directory where the results will be saved

%% Load images
% image1: The target image to which image2 will be registered.
% image2: The image to be transformed such that it is registered to image1.
image1_mean = double(read_tiff(session1_mean_path));
image2_mean = double(read_tiff(session2_mean_path));
image1_max = double(read_tiff(session1_max_path));
image2_max = double(read_tiff(session2_max_path));

%% Prepare inputs for PatchWarp
% Users can use any kinds of 2D summary images. Here, we are using a mean
% image and a max-projection image. You can use either only one of them, or
% use both by concatenating the images along the 3rd dimension. If you use
% both of them, the PatchWarp results are more robust and accurate. Users can
% add more summary images along the 3rd dimension if users need more robust 
% results(e.g. standard deviation map, correlation map).
image1_all = cat(3, image1_mean, image1_max);
image2_all = cat(3, image2_mean, image2_max);

%% Apply PatchWarp to register images from session#2 to images from session#1
% warp_blocksize:                       Row and column numbers for splitting FOV. Each image is split into [warp_blocksize]*[warp_blocksize] subfields 
%                                       for estimating and applying affine transformation matrices.
% warp_overlap_pix_frac:                Fraction of edge pixels that overlaps with the adjacent patches. [warp_overlap_pix_frac]*[length of a patch] pixels 
%                                       at the edge of each patch will be shared with the adjacent patch.
%                                       with its adjacent subfields.
% warp_init_1st_affine:                 Initial guess of affine transformation matrix for the 1st iteration of affine transformation which uses the who FOV. 
%                                       Default is the identity matrix (Note that the 3rd row is omitted). If the result is not satisfactory and you can guess 
%                                       the transfomation matrix (translation, rotation, etc.), changing this initial guess of the transformation matrix may improve the result.
warp_blocksize = 8;
warp_overlap_pix_frac = 0.15;
warp_init_1st_affine = [1 0 0; 0 1 0];
patchwarp_results = patchwarp_across_sessions(image1_all, image2_all, warp_blocksize, warp_overlap_pix_frac, warp_init_1st_affine);

save(fullfile(save_path, 'patchwarp_across_session_results.mat'), 'patchwarp_results')

%% Plot registered images
figure
set(gcf,'Position',[50 50 1550 950])
subplot(2,3,1)
imshow(imfuse(patchwarp_results.image1_all(:, :, 1),patchwarp_results.image2_all(:, :, 1),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
title('Mean images (raw)')
subplot(2,3,2)
imshow(imfuse(patchwarp_results.image1_all(:, :, 1),patchwarp_results.image2_affine1(:, :, 1),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
title('Mean images after 1st affine')
subplot(2,3,3)
imshow(imfuse(patchwarp_results.image1_all(:, :, 1),patchwarp_results.image2_affine2(:, :, 1),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
title('Mean images after 2nd affine')
subplot(2,3,4)
imshow(imfuse(patchwarp_results.image1_all(:, :, 2),patchwarp_results.image2_all(:, :, 2),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
title('Max images (raw)')
subplot(2,3,5)
imshow(imfuse(patchwarp_results.image1_all(:, :, 2),patchwarp_results.image2_affine1(:, :, 2),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
title('Max images after 1st affine')
subplot(2,3,6)
imshow(imfuse(patchwarp_results.image1_all(:, :, 2),patchwarp_results.image2_affine2(:, :, 2),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
title('Max images after 2nd affine')

%% Apply the obtained transformations to other images (e.g. individual frames, ROI mask)
% INPUTS:
% patchwarp_results:    Output from the above "patchwarp_across_sessions" function
% input_images:         Input can be 2D (x, y) or 3D (x, y, frames). For example, it can transform individual frames or a single summary image (e.g. ROI mask).    
% second_affine:        Apply only the 1st whole-FOV affine only (0), or apply both the 1st and the 2nd subfield-wise affine transformations (1).
% OUTPUTS:
% images_affine1:    	Images after the 1st affine transformation.
% images_affine2:       Images after the 2nd affine transformation. This will be empty if second_affine==0.

% image_path = 'G:\Data\171026\RH825\post_warp_affine\downsampled\roi_mask.tif';
image_path = 'G:\Data\171026\RH825\post_warp_affine\RH825camk2GC6sx18p70MatchingRSCr_001_002_corrected_warped.tif';
second_affine = 1;

input_images = double(read_tiff(image_path));
[images_affine1, images_affine2] = patchwarp_across_sessions_apply(input_images, patchwarp_results, second_affine);
figure
set(gcf,'Position',[50 50 1550 950])
subplot(2,3,1)
imshow(imfuse(input_images(:, :, 1), images_affine1(:, :, 1),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
