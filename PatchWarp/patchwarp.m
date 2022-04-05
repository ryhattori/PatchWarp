function patchwarp(ops)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PatchWarp
% -------------------
% 
% Released by Ryoma Hattori
% Email: rhattori0204@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if all input parameters exist
if ~isfield(ops, 'source_path')
    error('Specify ops.source_path!')
end
if ~isfield(ops, 'save_path')
    error('Specify ops.save_path!')
end
if ~isfield(ops, 'n_ch')
    disp('ops.n_ch is not specified. Running with ops.n_ch = 1')
    ops.n_ch = 1;
end
if ~isfield(ops, 'align_ch')
    disp('ops.align_ch is not specified. Running with ops.align_ch = 1')
    ops.align_ch = 1;
end
if ~isfield(ops, 'save_ch')
    disp('ops.save_ch is not specified. Running with ops.save_ch = 1')
    ops.save_ch = 1;
end
if ~isfield(ops, 'run_rigid_mc')
    disp('ops.run_rigid_mc is not specified. Running with ops.run_rigid_mc = 1')
    ops.run_rigid_mc = 1;
end
if ~isfield(ops, 'run_affine_wc')
    disp('ops.run_affine_wc is not specified. Running with ops.run_affine_wc = 1')
    ops.run_affine_wc = 1;
end
if ~isfield(ops, 'rigid_norm_method')
    disp('ops.rigid_norm_method is not specified. Running with ops.rigid_norm_method = rank')
    ops.rigid_norm_method = 'rank';
end
if ~isfield(ops, 'rigid_norm_radius')
    ops.rigid_norm_radius = '32';
end
if ~isfield(ops, 'rigid_template_block_num')
    disp('ops.rigid_template_block_num is not specified. Running with ops.rigid_template_block_num = 3')
    ops.rigid_template_block_num = 3;
end
if ~isfield(ops, 'rigid_template_threshold')
    disp('ops.rigid_template_threshold is not specified. Running with ops.rigid_template_threshold = 0.2')
    ops.rigid_template_threshold = 0.2;
end
if ~isfield(ops, 'rigid_template_tiffstack_num')
    disp('ops.rigid_template_tiffstack_num is not specified. Running with ops.rigid_template_tiffstack_num = 3')
    ops.rigid_template_tiffstack_num = 3;
end
if ~isfield(ops, 'rigid_template_center_frac')
    disp('ops.rigid_template_center_frac is not specified. Running with ops.rigid_template_center_frac = 0.8')
    ops.rigid_template_center_frac = 0.8;
end
if ~isfield(ops, 'affine_norm_radius')
    disp('ops.affine_norm_radius is not specified. Running with ops.affine_norm_radius = 32')
    ops.affine_norm_radius = 32;
end
if ~isfield(ops, 'warp_template_tiffstack_num')
    disp('ops.warp_template_tiffstack_num is not specified. Running with ops.warp_template_tiffstack_num = 5')
    ops.warp_template_tiffstack_num = 5;
end
if ~isfield(ops, 'warp_movave_tiffstack_num')
    disp('ops.warp_movave_tiffstack_num is not specified. Running with ops.warp_movave_tiffstack_num = 1')
    ops.warp_movave_tiffstack_num = 1;
end
if ~isfield(ops, 'warp_blocksize')
    disp('ops.warp_blocksize is not specified. Running with ops.warp_blocksize = 8')
    ops.warp_blocksize = 8;
end
if ~isfield(ops, 'warp_overlap_pix_frac')
    disp('ops.warp_overlap_pix_frac is not specified. Running with ops.warp_overlap_pix_frac = 0.15')
    ops.warp_overlap_pix_frac = 0.15;
end
if ~isfield(ops, 'n_split4warpinit')
    disp('ops.n_split4warpinit is not specified. Running with ops.n_split4warpinit = 6')
    ops.n_split4warpinit = 6;
end
if ~isfield(ops, 'edge_remove_pix')
    ops.edge_remove_pix = 0;
end
if ~isfield(ops, 'affinematrix_abssum_threshold')
    ops.affinematrix_abssum_threshold = 50;
end
if ~isfield(ops, 'affinematrix_abssum_jump_threshold')
    ops.affinematrix_abssum_jump_threshold = 10;
end
if ~isfield(ops, 'affinematrix_rho_threshold')
    ops.affinematrix_rho_threshold = 0.5;
end
if ~isfield(ops, 'affinematrix_medfilt_tiffstack_num')
    disp('ops.affinematrix_medfilt_tiffstack_num is not specified. Running with ops.affinematrix_medfilt_tiffstack_num = 1 (no filtering). Use larger value if warp result is unstable.')
    ops.affinematrix_medfilt_tiffstack_num = 1;
end
if ~isfield(ops, 'transform')
    disp('ops.transform is not specified. Running with ops.transform = affine')
    ops.transform = 'affine';
end
if ~isfield(ops, 'warp_pyramid_levels')
    ops.warp_pyramid_levels = 1;
end
if ~isfield(ops, 'warp_pyramid_iterations')
    disp('ops.warp_pyramid_iterations is not specified. Running with ops.warp_pyramid_iterations = 50')
    ops.warp_pyramid_iterations = 50;
end
if ~isfield(ops, 'downsample_frame_num')
    disp('ops.downsample_frame_num is not specified. Running with ops.downsample_frame_num = 50')
    ops.downsample_frame_num = 50;
end
if ~isfield(ops, 'worker_num')
    disp('ops.worker_num is not specified. Running with maximum worker number. Limit the number if you see out of memory error.')
    ops.worker_num = parcluster('local').NumWorkers;
end
if ~isfield(ops, 'network_temp_copy')
    ops.network_temp_copy = 1;
end

%% Start parallel pool
delete(gcp('nocreate'));
parpool(ops.worker_num);

%% Run rigid correction
save_path_prewarp = fullfile(ops.save_path, 'pre_warp');
save_path_postwarp = fullfile(ops.save_path, 'post_warp');
if ops.run_rigid_mc == 1
    disp('Performing rigid motion correction...')
    patchwarp_rigid(ops.source_path, save_path_prewarp, ops.n_ch, ops.align_ch, ops.save_ch, ops.rigid_norm_method, ops.rigid_norm_radius, ops.rigid_template_block_num, ops.rigid_template_threshold, ops.rigid_template_tiffstack_num, ops.rigid_template_center_frac, ops.downsample_frame_num, ops.network_temp_copy, []);
else
    disp('Skipping rigid motion correction...')
end

%% Run warp correction
if ops.run_affine_wc == 1
    disp('Performing warp correction...')
    patchwarp_affine(save_path_prewarp, save_path_postwarp, ops.n_ch, ops.align_ch, ops.affine_norm_radius, ops.warp_template_tiffstack_num, ops.warp_movave_tiffstack_num, ops.warp_blocksize, ops.warp_overlap_pix_frac, ops.n_split4warpinit, ops.edge_remove_pix, ops.affinematrix_abssum_threshold, ops.affinematrix_abssum_jump_threshold, ops.affinematrix_rho_threshold, ops.affinematrix_medfilt_tiffstack_num, ops.transform, ops.warp_pyramid_levels, ops.warp_pyramid_iterations, ops.downsample_frame_num, ops.worker_num, ops.network_temp_copy);
else
    disp('Skipping warp correction...')
end

%% Shutdown parallel pool
delete(gcp('nocreate'))
        
    
    
    
    