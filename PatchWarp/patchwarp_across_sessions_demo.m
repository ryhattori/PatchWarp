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
session1_mean_path = '...\early_session_mean_image.tif';
session1_max_path = '...\early_session_max_projection_image.tif';
session2_mean_path = '...\late_session_mean_image.tif';
session2_max_path = '...\late_session_max_projection_image.tif';

%% Specify saving directory
save_path = 'G:\Data\171026\RH825';    % Directory where the results will be saved

%% Load images
% image1: The target image to which image2 will be registered.
% image2: The image to be transformed such that it is registered to image1.
image1_mean = double(read_tiff(session1_mean_path));
image2_mean = double(read_tiff(session2_mean_path));
image1_max = double(read_tiff(session1_max_path));
image2_max = double(read_tiff(session2_max_path));

%% Normalize image intensity
image1_mean = 1000*reshape(normalize(image1_mean(:), 'range'), size(image1_mean));
image2_mean = 1000*reshape(normalize(image2_mean(:), 'range'), size(image2_mean));
image1_max = 1000*reshape(normalize(image1_max(:), 'range'), size(image1_max));
image2_max = 1000*reshape(normalize(image2_max(:), 'range'), size(image2_max));
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
% transform1:                           Transformation type for the 1st transformation that uses the whole FOV. Typically 'euclidean' should work best. 
%                                       If 'euclidean' does not work, try 'affine'. This 1st transformation matrix is obtained by the pyramid method with 3 levels.
% transform2:                           Transformation type for the 2nd transformations for each subfield. Typically 'affine' should work best. 
%                                       These 2nd transformation matrices are obtained without pyramids.
% warp_blocksize:                       Row and column numbers for splitting FOV. Each image is split into [warp_blocksize]*[warp_blocksize] subfields 
%                                       for estimating and applying affine transformation matrices. 
% warp_overlap_pix_frac:                Fraction of edge pixels that overlaps with the adjacent patches. [warp_overlap_pix_frac]*[length of a patch] pixels 
%                                       at the edge of each patch will be shared with the adjacent patch.
%                                       with its adjacent subfields.
% norm_radius:                          If this is > 0, the intensity of each pixel is normalized with the mean intensity within the specified radius. 
%                                       Disable this normalization by setting this to 0. Disable when the signals are sparse (e.g. dendrite imaging, axon imaging).
%                                       
transform1 = 'euclidean';
transform2 = 'affine';
warp_blocksize = 6;
warp_overlap_pix_frac = 0.15;
norm_radius = 32;     % Set to 0 when the signals are sparse or the result does not look good.
patchwarp_results = patchwarp_across_sessions(image1_all, image2_all, transform1, transform2, warp_blocksize, warp_overlap_pix_frac, norm_radius);

save(fullfile(save_path, 'patchwarp_across_session_results.mat'), 'patchwarp_results')

%% Plot registered images
figure
set(gcf,'Position',[50 50 1550 950])
subplot(2,3,1)
imshow(imfuse(patchwarp_results.image1_all(:, :, 1),patchwarp_results.image2_all(:, :, 1),'falsecolor','Scaling','joint','ColorChannels','red-cyan'));
title('Mean images (raw)')
subplot(2,3,2)
imshow(imfuse(patchwarp_results.image1_all(:, :, 1),patchwarp_results.image2_warp1(:, :, 1),'falsecolor','Scaling','joint','ColorChannels','red-cyan'));
title('Mean images after 1st transformation (euclidean)')
subplot(2,3,3)
imshow(imfuse(patchwarp_results.image1_all(:, :, 1),patchwarp_results.image2_warp2(:, :, 1),'falsecolor','Scaling','joint','ColorChannels','red-cyan'));
title('Mean images after 2nd transformation (affine)')
subplot(2,3,4)
imshow(imfuse(patchwarp_results.image1_all(:, :, 2),patchwarp_results.image2_all(:, :, 2),'falsecolor','Scaling','joint','ColorChannels','red-cyan'));
title('Max images (raw)')
subplot(2,3,5)
imshow(imfuse(patchwarp_results.image1_all(:, :, 2),patchwarp_results.image2_warp1(:, :, 2),'falsecolor','Scaling','joint','ColorChannels','red-cyan'));
title('Max images after 1st transformation (euclidean)')
subplot(2,3,6)
imshow(imfuse(patchwarp_results.image1_all(:, :, 2),patchwarp_results.image2_warp2(:, :, 2),'falsecolor','Scaling','joint','ColorChannels','red-cyan'));
title('Max images after 2nd transformation (affine)')

%% Apply the obtained transformations to other images (e.g. individual frames, ROI mask)
% INPUTS:
% patchwarp_results:    Output from the above "patchwarp_across_sessions" function
% input_images:         Input can be 2D (x, y) or 3D (x, y, frames). For example, it can transform individual frames or a single summary image (e.g. ROI mask).    
% second_affine:        Apply only the 1st whole-FOV euclidean only (0), or apply both the 1st and the 2nd subfield-wise affine transformations (1).
% OUTPUTS:
% images_warp1:         Images after the 1st transformation (euclidean).
% images_warp2:         Images after the 2nd transformation (affine). This will be empty if second_affine==0.

second_affine = 1;
image_path = 'G:\Data\171026\RH825\corrected\post_warp_affine\downsampled\roi_mask.tif';
% image_path = 'G:\Data\171026\RH825\corrected\post_warp_affine\RH825camk2GC6sx18p70MatchingRSCr_001_002_corrected_warped.tif';

input_images = double(read_tiff(image_path));

[images_warp1, images_warp2] = patchwarp_across_sessions_apply(input_images, patchwarp_results, second_affine);

% An input image can be any images with the same size from session#2. e.g. Binary cellular ROI mask, 3D tiff stack, summary image.
% For example, if '...\late_session_max_projection_image.tif'is used as the input_image, it will return the same transformed image as in
% "patchwarp_results.image2_warp1/2" from the previous section. Note that the edge pixels have 0 values in these warped images.

