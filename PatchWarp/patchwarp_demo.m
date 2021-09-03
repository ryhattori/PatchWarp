%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PatchWarp demo
% -------------------
% 
% Released by Ryoma Hattori
% Email: rhattori0204@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add PatchWarp directory to MATLAB path
patchwarp_path = 'Z:\People\Ryoma\PatchWarp';
addpath(genpath(patchwarp_path))

%% Specifiy source directory and saving directory
source_path = 'Z:\People\Ryoma\MC_data\171104\RH825';    % Directory that contains original tif stack files
save_path = 'Z:\People\Ryoma\MC_data\171104\RH825\corrected';    % Directory where motion corrected images will be saved

%% Set general parameters
% n_ch:                     Number of saved PMT channels. Set this to 2 for if the session was 2-color imaging.
% align_ch:                 Imaging channel to be aligned
% save_ch:                  Imaging channel to be saved
% run_rigid_mc:             skip rigid motion correction (0), or do ridig motion correction (1)
% run_affine_wc:            skip warp correction (0), or do warp correction(1)
n_ch = 1;
align_ch = 1;
save_ch = 1;
run_rigid_mc = 1;
run_affine_wc = 1;

%% Set parameter for rigid motion correction
% rigid_template_tiffstack_num: Number of tif stack files used for estimating template images. For example, if each tif stack file has 500
%                               frames, each template image will be created using 500*[rigid_template_tiffstack_num] frames.
% rigid_template_block_num:     Number of blocks for an imaging session to be split for rigid motion correction. 
%                               Template image for rigid motion correction is re-estimated for each block.
% rigid_template_threshold:     Quantile threshold for selecting the frames that are used for template images. Pearson correlation coefficient is
%                               calculated between the mean of all frames within a window (specified by rigid_template_tiffstack_num) and each 
%                               frame within the window. Frames with the correlation coefficient lower than the specified quantile threshold will
%                               be ignored when obtaining the final template image.
% rigid_template_center_frac:   Central fraction of a template image used for registration. For example, 5% of pixels are ignored from each edge
%                               when rigid_template_center_frac = 0.9.
rigid_template_tiffstack_num = 5;
rigid_template_block_num = 5;   % This must be an odd number. Minimum is 3.
rigid_template_threshold = 0.2;
rigid_template_center_frac = 0.95;

%% Set parameter for warp correction
% transform:                            Type of the image transformation. 'affine' works best. ('translation', 'euclidean', 'affine', 'homography')
% warp_pyramid_levels:                  The number of levels in pyramid scheme (1 for a non pyramid implementation).
% warp_pyramid_iterations:              The number of iterations per level.
% warp_template_tiffstack_num:          Number of tif stack files that are used for making the template image.
% warp_movave_tiffstack_num:            Window size for temporally smoothing downsampled images before estimation of affine transformation matrices.
% warp_blocksize:                       Row and column numbers for splitting FOV. Each image is split into [warp_blocksize]*[warp_blocksize] subfields 
%                                       for estimating and applying affine transformation matrices.
% warp_overlap_pix_frac:                Fraction of edge pixels that overlaps with the adjacent patches. [warp_overlap_pix_frac]*[length of a patch] pixels 
%                                       at the edge of each patch will be shared with the adjacent patch.
%                                       with its adjacent subfields.
% edge_remove_pix:                      Number of pixels that you want to ignore from the lateral edges of a FOV. This may be useful when your FOV contains
%                                       edges that were not scanned by a microscope. 
% n_split4warpinit:                     Sessions are split to this number, and initial affine transformation matrices are re-estimated for each block.
%                                       Gradient descent algorithm find the best matrices from the initial guess.
% affinematrix_abssum_threshold:        If the sum of absolute values of all elements in an estimated affine transformation matrix exceeds 
%                                       affinematrix_abssum_threshold, the matrix will be ignored before median temporal filtering.
% affinematrix_abssum_jump_threshold:   Ignore estimated affine matrices when the estimated matrix is very different from the matrix from the adjacent frame. 
%                                       If the sum(abs(matrix1 - matrix2)) > affinematrix_abssum_jump_threshold, the matrix will be ignored.
% affinematrix_rho_threshold:           If the enhanced correlation coefficient between the template and transformed subfield is less than or equal to
%                                       this threshold, the matrix will be ignored before median temporal filtering.
% affinematrix_medfilt_tiffstack_num:   Window size (number of tif stack files) that are used for median temporal filtering of affine transformation matrices.
transform = 'affine';
warp_pyramid_levels = 1;
warp_pyramid_iterations = 50;
warp_template_tiffstack_num = 11;
warp_movave_tiffstack_num = 21;
warp_blocksize = 8;     % For moderate distortion, use small number (e.g. 2-4). For severe distortion, use large number. Note that the processing time takes
                        % much longer if you use a large blocksize.
warp_overlap_pix_frac = 0.15;
edge_remove_pix = 0;
n_split4warpinit = 6;   % This must be an even number.
affinematrix_abssum_threshold = 50;
affinematrix_abssum_jump_threshold = 10;
affinematrix_rho_threshold = 0.5;
affinematrix_medfilt_tiffstack_num = 15;

%% Set parameter for a downsampled motion corrected tiff stack
% This tiff stack file can be used for visual inspection of the motion correction quality
% downsample_frame_num:     Window size for non-overlapping moving averaging.
downsample_frame_num = 50; 

%% Run PatchWarp
patchwarp(source_path, save_path, n_ch, align_ch, save_ch, run_rigid_mc, run_affine_wc,...
    rigid_template_block_num, rigid_template_threshold, rigid_template_tiffstack_num, rigid_template_center_frac,...
    warp_template_tiffstack_num, warp_movave_tiffstack_num, warp_blocksize, warp_overlap_pix_frac, n_split4warpinit, edge_remove_pix,...
    affinematrix_abssum_threshold, affinematrix_abssum_jump_threshold, affinematrix_rho_threshold, affinematrix_medfilt_tiffstack_num, transform, warp_pyramid_levels, warp_pyramid_iterations,...
    downsample_frame_num);


