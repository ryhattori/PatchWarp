function patchwarp_rigid(source_path, save_path, n_ch, align_ch, save_ch, rigid_norm_radius, rigid_template_block_num, rigid_template_threshold, rigid_template_tiffstack_num, rigid_template_center_frac, n_downsampled, n_downsampled_perstack, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PatchWarp
% patchwarp_rigid
% -------------------
% 
% Iteratively re-estimates template images for rigid motion correction
% The core rigid motion correction algorithm (Pyramid approach) was written by Aki Mitani and is available from https://github.com/amitani/matlab_motion_correct
% Details of the Pyramid approach is described in 
% Mitani & Komiyama, Real-Time Processing of Two-Photon Calcium Imaging Data Including Lateral Motion Artifact Correction., Front. Neuroinform., 18, 2018
%
% Released by Ryoma Hattori
% Email: rhattori0204@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(~exist('opt','var'))
        opt='';
    end
    if(~exist('align_ch','var'))
        align_ch=[];
    end
    if(~exist('save_ch','var'))
        save_ch=[];
    end
    
    L=Logger;
    
    L.newline('Initializing');
    
    if(~exists_file(save_path))
        mkdir(save_path);
    end
    
    [~, fn_list] = fastdir(source_path, common_regexp('tiff_ext'));

    if isempty(n_downsampled_perstack)
        n_downsampled_perstack = inf;
    end
    
    middle_tiffstackfile_id = round(length(fn_list)/2);
    target_save_path = fullfile(save_path,'target');
    target_fn = cell(rigid_template_block_num, 1);
    for i = 1:rigid_template_block_num 
    	target_fn{i} = fullfile(target_save_path, ['template_AVG', num2str(i), '.tif']);
    end
    if(~exists_file(target_fn{1}))
        L.newline('Making target from %s', fn_list{middle_tiffstackfile_id});
        target =  make_template_from_file_multiple(fn_list, middle_tiffstackfile_id - floor((rigid_template_tiffstack_num - 1)/2):middle_tiffstackfile_id + ceil((rigid_template_tiffstack_num - 1)/2),align_ch, n_ch, rigid_template_threshold, false); % 4th is the threshold for quantile(c,1-threshold). 1: use all frames to make a template. 
        target = parse_image_input(target,align_ch);
        if(~exists_file(save_path))
            L.newline('Make save dir');
            mkdir(save_path);
        end
        if(~exists_file(target_save_path))
            mkdir(target_save_path);
        end
        
        parfor i = middle_tiffstackfile_id - floor((rigid_template_tiffstack_num - 1)/2):middle_tiffstackfile_id + ceil((rigid_template_tiffstack_num - 1)/2)
            pyramid_registration(fn_list{i}, target, target_save_path, align_ch, save_ch, n_downsampled, n_downsampled_perstack, n_ch, rigid_norm_radius, rigid_template_center_frac);      
        end
        fn_list_corrected_temp1 = cell(rigid_template_tiffstack_num, 1);
        for i = 1:rigid_template_tiffstack_num
            [~,fn_list_corrected_temp1{i},~] = fileparts(fn_list{i - 1 + middle_tiffstackfile_id - floor((rigid_template_tiffstack_num - 1)/2)});
            fn_list_corrected_temp1{i} = fullfile(target_save_path,[fn_list_corrected_temp1{i} '_corrected.tif']);
        end
        target = make_template_from_file_multiple(fn_list_corrected_temp1, 1:rigid_template_tiffstack_num, 1, 1, rigid_template_threshold, false);
        write_tiff(target_fn{1}, int16(target));
    end
    processed = false;
    
    fn_list_corrected = cell(size(fn_list,1),1);
    for i = 1:size(fn_list)
        [~,fn_list_corrected{i},~] = fileparts(fn_list{i});
        fn_list_corrected{i} = fullfile(save_path,[fn_list_corrected{i} '_corrected.tif']);
    end
    
    if((strcmp(opt,'f')||strcmp(opt,'force')||~pyramid_registration(fn_list{end},[], save_path, align_ch, save_ch, n_downsampled, n_downsampled_perstack, n_ch))) && ~strcmp(opt,'a')
        block_range_list = zeros(rigid_template_block_num, 2);
        block_range_list(1, 1) = floor(rigid_template_block_num/2) * round(length(fn_list)/rigid_template_block_num);
        block_range_list(1, 2) = ceil(rigid_template_block_num/2) * round(length(fn_list)/rigid_template_block_num);
        block_range_list(rigid_template_block_num - 1, 1) = (rigid_template_block_num - 1) * round(length(fn_list)/rigid_template_block_num) + 1;
        block_range_list(rigid_template_block_num - 1, 2) = length(fn_list);
        block_range_list(rigid_template_block_num, 1) = 1;
        block_range_list(rigid_template_block_num, 2) = round(length(fn_list)/rigid_template_block_num) - 1;
        if rigid_template_block_num > 3
            if ~mod(rigid_template_block_num,2)
               error('rigid_template_block_num should be an odd number!!!') 
            end
            for i = 1:(rigid_template_block_num - 3)/2
                block_range_list(2 * i, 1) = (ceil(rigid_template_block_num/2) + i - 1) * round(length(fn_list)/rigid_template_block_num) + 1;
                block_range_list(2 * i, 2) = (ceil(rigid_template_block_num/2) + i) * round(length(fn_list)/rigid_template_block_num);
                block_range_list(2 * i + 1, 1) = (floor(rigid_template_block_num/2) - i) * round(length(fn_list)/rigid_template_block_num);
                block_range_list(2 * i + 1, 2) = (floor(rigid_template_block_num/2) - i + 1) * round(length(fn_list)/rigid_template_block_num) - 1;
            end
        end
        
        block_id = 1;
        parfor stack_id = block_range_list(block_id, 1):block_range_list(block_id, 2)
            pyramid_registration(fn_list{stack_id}, target_fn{block_id}, save_path, align_ch, save_ch, n_downsampled, n_downsampled_perstack, n_ch, rigid_norm_radius, rigid_template_center_frac);      
        end
        target = make_template_from_file_multiple(fn_list_corrected, block_range_list(block_id, 2) - rigid_template_tiffstack_num + 1:block_range_list(block_id, 2), align_ch, n_ch, rigid_template_threshold, false);
        target = parse_image_input(target, align_ch);
        write_tiff(target_fn{2}, int16(target));
        target = make_template_from_file_multiple(fn_list_corrected, block_range_list(block_id, 1):block_range_list(block_id, 1) + rigid_template_tiffstack_num - 1, align_ch, n_ch, rigid_template_threshold, false);
        target = parse_image_input(target, align_ch);
        write_tiff(target_fn{3}, int16(target));
        
        for i = 1:(rigid_template_block_num - 1)/2
            parfor stack_id = block_range_list(2 * i, 1):block_range_list(2 * i, 2)
                pyramid_registration(fn_list{stack_id}, target_fn{2 * i}, save_path, align_ch, save_ch, n_downsampled, n_downsampled_perstack, n_ch, rigid_norm_radius, rigid_template_center_frac);      
            end
            if i ~=(rigid_template_block_num - 1)/2
                target = make_template_from_file_multiple(fn_list_corrected, block_range_list(2 * i, 2) - rigid_template_tiffstack_num + 1:block_range_list(2 * i, 2), align_ch, n_ch, rigid_template_threshold, false);
                target = parse_image_input(target, align_ch);
                write_tiff(target_fn{2 * (i + 1)}, int16(target));
            end
            parfor stack_id = block_range_list(2 * i + 1, 1):block_range_list(2 * i + 1, 2)
                pyramid_registration(fn_list{stack_id}, target_fn{2 * i + 1}, save_path, align_ch, save_ch, n_downsampled, n_downsampled_perstack, n_ch, rigid_norm_radius, rigid_template_center_frac);      
            end
            if i ~= (rigid_template_block_num - 1)/2
                target = make_template_from_file_multiple(fn_list_corrected, block_range_list(2 * i + 1, 1):block_range_list(2 * i + 1, 1) + rigid_template_tiffstack_num - 1, align_ch, n_ch, rigid_template_threshold, false);
                target = parse_image_input(target, align_ch);
                write_tiff(target_fn{2 * (i + 1) + 1}, int16(target));
            end
        end
            
        processed = true;
    else
        processed = true;
        disp('skipping because already motion corrected.');
    end
    if(processed || ~make_downsampled_tiff(save_path,'-s') )
        make_downsampled_tiff(save_path);
    end
    
    L.newline('done');
end

function b = exists_file(fn)
    b = java.io.File(fn).exists()>0;
end
