function patchwarp(source_path, save_path, n_ch, align_ch, save_ch, run_rigid_mc, run_affine_wc, rigid_norm_radius, rigid_template_block_num, rigid_template_threshold, rigid_template_tiffstack_num, rigid_template_center_frac, affine_norm_radius, warp_template_tiffstack_num, warp_movave_tiffstack_num, warp_blocksize, warp_overlap_pix_frac, n_split4warpinit, edge_remove_pix, affinematrix_abssum_threshold, affinematrix_abssum_jump_threshold, affinematrix_rho_threshold, affinematrix_medfilt_tiffstack_num, transform, warp_pyramid_levels, warp_pyramid_iterations, downsample_frame_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PatchWarp
% -------------------
% 
% Released by Ryoma Hattori
% Email: rhattori0204@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_path_prewarp = fullfile(save_path, 'pre_warp');
save_path_postwarp = fullfile(save_path, 'post_warp');
if run_rigid_mc == 1
    disp('Performing rigid motion correction...')
    patchwarp_rigid(source_path, save_path_prewarp, n_ch, align_ch, save_ch, rigid_norm_radius, rigid_template_block_num, rigid_template_threshold, rigid_template_tiffstack_num, rigid_template_center_frac, downsample_frame_num, []);
else
    disp('Skipping rigid motion correction...')
end

if run_affine_wc == 1
    disp('Performing warp correction...')
    patchwarp_affine(save_path_prewarp, save_path_postwarp, n_ch, align_ch, affine_norm_radius, warp_template_tiffstack_num, warp_movave_tiffstack_num, warp_blocksize, warp_overlap_pix_frac, n_split4warpinit, edge_remove_pix, affinematrix_abssum_threshold, affinematrix_abssum_jump_threshold, affinematrix_rho_threshold, affinematrix_medfilt_tiffstack_num, transform, warp_pyramid_levels, warp_pyramid_iterations, downsample_frame_num);
else
    disp('Skipping warp correction...')
end
        
    
    
    
    