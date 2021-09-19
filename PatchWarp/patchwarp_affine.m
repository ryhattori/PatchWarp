function patchwarp_affine(source_path, save_path, n_ch, align_ch, affine_norm_radius, warp_template_tiffstack_num, warp_movave_tiffstack_num, warp_blocksize, warp_overlap_pix_frac, n_split4warpinit, edge_remove_pix, affinematrix_abssum_threshold, affinematrix_abssum_jump_threshold, affinematrix_rho_threshold, affinematrix_medfilt_tiffstack_num, transform, warp_pyramid_levels, warp_pyramid_iterations, downsample_frame_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PatchWarp
% pacman_warp
% -------------------
% Piece-wise affine transformation of rigid-motion corrected images
% 
% Released by Ryoma Hattori
% Email: rhattori0204@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if ispc
%    %make temp folder and files
%    local_dir = 'C:\temp_files';
%    local_filename = [local_dir filesep 'temp_img_' datestr(now,'yymmddHHMMSSFFF') '.tif'];
%    if ~isdir(local_dir)
%        mkdir(local_dir)
%    end
%    % delete local temp file if it exists
%    if exist(local_filename,'file')
%        delete(local_filename)
%    end
% end

mkdir(save_path, 'downsampled');

%find downsampled_perstack file % get names to save
fn_downsampled_perstack = dir(fullfile(source_path,'downsampled','downsampled_perstack.tif'));
fn_downsampled = dir(fullfile(source_path,'downsampled',['*downsampled_', num2str(downsample_frame_num), '.tif']));
if length(fn_downsampled_perstack)~=1 || length(fn_downsampled)~=1
    disp(fn_downsampled_perstack); disp(fn_downsampled);
    error('downsampled_perstack file is not the only one')
end
downsampled_perstack_name = fn_downsampled_perstack.name;
fn_downsampled_name = fn_downsampled.name;
fn_downsampled_perstack = fullfile(source_path,'downsampled',downsampled_perstack_name);
fn_downsampled_perstack_save = fullfile(save_path,'downsampled',[downsampled_perstack_name(1:end-4),'_warped.tif']);
fn_affinetransmat = fullfile(save_path, 'affine_transformation_matrix.mat');
% fn_downsampled = fullfile(source_path,'downsampled',fn_downsampled_name);
fn_downsampled2save = fullfile(save_path,'downsampled',[fn_downsampled_name(1:end-4),'_warped.tif']);
fn_downsampled2save_max = fullfile(save_path,'downsampled',[fn_downsampled_name(1:end-4),'_warped_max.tif']);
fn_downsampled2save_mean = fullfile(save_path,'downsampled',[fn_downsampled_name(1:end-4),'_warped_mean.tif']);
% fn_downsampled2save_max_norm = fullfile(save_path,'downsampled',[fn_downsampled_name(1:end-4),'_warped_max_normalized.tif']);
% fn_downsampled2save_mean_norm = fullfile(save_path,'downsampled',[fn_downsampled_name(1:end-4),'_warped_mean_normalized.tif']);

%default parameters for warping
% n_ch = 1;
% ch = 1;
% margin_ratio = 0.125;
% warp_pyramid_levels = 1;
% warp_pyramid_iterations = 50;
% transform = 'affine';

disp('Loading downsampled_perstack file...')
[stack_downsampled_perstack, info_downsampled_perstack] = read_tiff(fn_downsampled_perstack, align_ch, n_ch);
% if ispc %copy to local
%     disp('Copying downsampled_perstack to local drive...')
%     img_filename = fn_downsampled_perstack;
%     disp(img_filename);
%     copyfile(img_filename,local_filename);
%     disp('Done.')
%     img_filename = local_filename;
%     disp(img_filename);
%     [stack_downsampled_perstack, info_downsampled_perstack] = read_tiff(img_filename, align_ch, n_ch);
% else
%     [stack_downsampled_perstack, info_downsampled_perstack] = read_tiff(fn_downsampled_perstack, align_ch, n_ch);
% end

stack_downsampled_perstack = stack_downsampled_perstack(:, 1 + edge_remove_pix:end - edge_remove_pix, :); % remove n pixels from left and right edges. 

zero_zone = min(stack_downsampled_perstack(:, :, 2:end-1),[],3); % zero-zone created by motion correction, but excluding 1st and last
nonzero_row = sum(zero_zone,2)~=0;
nonzero_column = sum(zero_zone,1)~=0;
stack_downsampled_perstack = stack_downsampled_perstack(nonzero_row, nonzero_column, :); % remove zero-zone

% minpixvalue = min(stack_downsampled_perstack(:,:,2:end-1),[],3);
% prctile(minpixvalue(:),50);

ny = size(stack_downsampled_perstack,1);
nx = size(stack_downsampled_perstack,2);
nz = size(stack_downsampled_perstack,3);

% warp_overlap_pix = round((overlap_pix-1)/2);
warp_overlap_pix = ceil(warp_overlap_pix_frac * (nx/warp_blocksize));

%% x-y ranges of each patch
qN_x = cell(warp_blocksize, warp_blocksize);
for i1 = 1:warp_blocksize
    qN_x{i1,1} = 1:ceil(nx/warp_blocksize)+warp_overlap_pix;
    for i2 = 2:warp_blocksize-1
        qN_x{i1,i2} = (i2-1)*ceil(nx/warp_blocksize)-warp_overlap_pix:i2*ceil(nx/warp_blocksize)+warp_overlap_pix;
    end
    qN_x{i1,warp_blocksize} = (warp_blocksize-1)*ceil(nx/warp_blocksize)-warp_overlap_pix:nx;
end
qN_y = cell(warp_blocksize,warp_blocksize);
for i2 = 1:warp_blocksize
    qN_y{1,i2} = 1:ceil(ny/warp_blocksize)+warp_overlap_pix;
    for i1 = 2:warp_blocksize-1
        qN_y{i1,i2} = (i1-1)*ceil(ny/warp_blocksize)-warp_overlap_pix:i1*ceil(ny/warp_blocksize)+warp_overlap_pix;
    end
    qN_y{warp_blocksize,i2} = (warp_blocksize-1)*ceil(ny/warp_blocksize)-warp_overlap_pix:ny;
end

%% Get affine transformation matrices in warp_cell
disp('Estimating affine transformation matrices...')

warp_cell = cell(warp_blocksize, warp_blocksize, nz);

if mod(n_split4warpinit, 2) == 1
    error('Session split must be even number for warping')
end

template_im = mean(stack_downsampled_perstack(:, :, round(nz/2)-round((warp_template_tiffstack_num - 1)/2):round(nz/2)- round((warp_template_tiffstack_num - 1)/2) + warp_template_tiffstack_num - 1), 3); % default (:,:,1)
template_im_normalized = imnormalize2(template_im, affine_norm_radius);

smoothed_downsampled_perstack = NaN(ny, nx, nz);
for i = 1:nz
    if i <= round((warp_movave_tiffstack_num - 1) / 2)
        smoothed_downsampled_perstack(:,:,i) = int16(imnormalize2(mean(stack_downsampled_perstack(:, :, 2:warp_movave_tiffstack_num), 3), affine_norm_radius));
    elseif i > nz - round((warp_movave_tiffstack_num - 1) / 2)
        smoothed_downsampled_perstack(:,:,i) = int16(imnormalize2(mean(stack_downsampled_perstack(:, :, nz-warp_movave_tiffstack_num:nz-1), 3), affine_norm_radius));
    else
        smoothed_downsampled_perstack(:,:,i) = int16(imnormalize2(mean(stack_downsampled_perstack(:, :, i-round((warp_movave_tiffstack_num - 1) / 2):i-round((warp_movave_tiffstack_num - 1) / 2) + warp_movave_tiffstack_num - 1),3), affine_norm_radius));
    end
end

if strcmp(transform, 'translation')
    warp_init = [0; 0];
elseif strcmp(transform, 'homography')
    warp_init = [1 0 0; 0 1 0; 0 0 1];
else
    warp_init = [1 0 0; 0 1 0];
end
warp4template = cell(warp_blocksize, warp_blocksize, n_split4warpinit);
for i1 = 1:warp_blocksize
    for i2 = 1:warp_blocksize
        warp4template{i1,i2,n_split4warpinit/2} = warp_init;
        warp4template{i1,i2,1+n_split4warpinit/2} = warp_init;
    end
end
rho = NaN(warp_blocksize,warp_blocksize,nz);

warp_range = nan(n_split4warpinit,2);
for i = 1:n_split4warpinit
    if i ==1
        warp_range(i,:) = [1,i*ceil(nz/n_split4warpinit)];
    elseif i == n_split4warpinit
        warp_range(i,:) = [(i-1)*ceil(nz/n_split4warpinit)+1, nz];
    else
        warp_range(i,:) = [(i-1)*ceil(nz/n_split4warpinit)+1, i*ceil(nz/n_split4warpinit)];
    end
end

if isempty(gcp('nocreate'))
    parpool;
end
pctRunOnAll warning off

n_warp2use4template = 7;
for nth_warp = [n_split4warpinit/2:-1:1,n_split4warpinit/2+1:n_split4warpinit]
    for i3 = warp_range(nth_warp,1):warp_range(nth_warp,2)
        parfor i1 = 1:warp_blocksize
            for i2 = 1:warp_blocksize
                  [~, warp_cell{i1,i2,i3}, ~, rho(i1,i2,i3)]= ...
                    ecc_patchwarp(smoothed_downsampled_perstack(qN_y{i1,i2}, qN_x{i1,i2},i3),...
                    template_im_normalized(qN_y{i1,i2}, qN_x{i1,i2}),...
                    warp_pyramid_levels, warp_pyramid_iterations, transform, warp4template{i1, i2, nth_warp});
%                     warp_pyramid_levels, warp_pyramid_iterations, transform, warp4template{i1, i2, nth_warp}(1:2,1:3));
            end
        end
    end
    
    for i1 = 1:warp_blocksize
        for i2 = 1:warp_blocksize
            for i3 = warp_range(nth_warp,1):warp_range(nth_warp,2)
                if sum(sum(abs(warp_cell{i1,i2,i3}) > affinematrix_abssum_threshold))
                    rho(i1,i2,i3) = NaN;
                end
                if isempty(warp_cell{i1,i2,i3})
                    if strcmp(transform, 'translation')
                        warp_cell{i1,i2,i3} = NaN(2,1);
                    elseif strcmp(transform, 'homography')
                        warp_cell{i1,i2,i3} = NaN(3,3);
                    else
                        warp_cell{i1,i2,i3} = NaN(2,3);
                    end
                end
                if ~strcmp(transform, 'translation')
                    if warp_cell{i1,i2,i3}(1,1)<0.6 || warp_cell{i1,i2,i3}(2,2)<0.6 || warp_cell{i1,i2,i3}(1,1)>1.4 || warp_cell{i1,i2,i3}(2,2)>1.4 || rho(i1,i2,i3) <= affinematrix_rho_threshold
                        if strcmp(transform, 'homography')
                            warp_cell{i1,i2,i3} = NaN(3,3);
                        else
                            warp_cell{i1,i2,i3} = NaN(2,3);
                        end
                    end
                else
                    if rho(i1,i2,i3) <= affinematrix_rho_threshold
                        warp_cell{i1,i2,i3} = NaN(2,1);
                    end
                end
                if (size(warp_cell{i1,i2,i3},1) == 3) && strcmp(transform, 'affine')
                    warp_cell{i1,i2,i3}(3,:) = [];
                end
                if (size(warp_cell{i1,i2,i3},1) ~= 3) && strcmp(transform, 'homography')
                    warp_cell{i1,i2,i3} = NaN(3,3);
                end
            end
        end
    end
    
    if (nth_warp ~= 1 && nth_warp ~= n_split4warpinit)
        for i1 = 1:warp_blocksize
            for i2 = 1:warp_blocksize
                if nth_warp <= n_split4warpinit/2
                    temp_id = find(~isnan(sum(sum(cell2mat(warp_cell(i1,i2,warp_range(nth_warp,1):warp_range(nth_warp,2))),1),2))==1,n_warp2use4template,'first');
                    if ~isempty(temp_id)
                        warp4template{i1,i2,nth_warp-1} = nanmedian(cell2mat(warp_cell(i1,i2,temp_id + warp_range(nth_warp,1) - 1)),3);
                    else
                        warp4template{i1,i2,nth_warp-1} = warp4template{i1,i2,nth_warp};
                    end
                    if isempty(warp4template{i1,i2,nth_warp-1})
                        warp4template{i1,i2,nth_warp-1} = warp4template{i1,i2,nth_warp};
                    end
                else
                    temp_id = find(~isnan(sum(sum(cell2mat(warp_cell(i1,i2,warp_range(nth_warp,1):warp_range(nth_warp,2))),1),2))==1,n_warp2use4template,'last');
                    if ~isempty(temp_id)
                        warp4template{i1,i2,nth_warp+1} = nanmedian(cell2mat(warp_cell(i1,i2,temp_id + warp_range(nth_warp,1) - 1)),3);
                    else
                        warp4template{i1,i2,nth_warp+1} = warp4template{i1,i2,nth_warp};
                    end
                    if isempty(warp4template{i1,i2,nth_warp+1})
                        warp4template{i1,i2,nth_warp+1} = warp4template{i1,i2,nth_warp};
                    end
                end
            end
        end
    end
end

% Clean warp_cell with temporal median filtering
half_medfilt_length = round(affinematrix_medfilt_tiffstack_num / 2);
for i3 = 1:nz
    for i1 = 1:warp_blocksize
        for i2 = 1:warp_blocksize
            if i3 <= half_medfilt_length
                correlated_id = rho(i1,i2,1:2*half_medfilt_length) > affinematrix_rho_threshold;
                warp_cell{i1,i2,i3} = nanmedian(cell2mat(warp_cell(i1,i2,correlated_id)),3);
            elseif i3 >= nz-(half_medfilt_length-1)
                correlated_id = rho(i1,i2,end-(2*half_medfilt_length-1):end) > affinematrix_rho_threshold;
                warp_cell{i1,i2,i3} = nanmedian(cell2mat(warp_cell(i1,i2,(nz-(2*half_medfilt_length-1)-1)+find(correlated_id))),3);
            else
                correlated_id = rho(i1,i2,i3-half_medfilt_length:i3+half_medfilt_length) > affinematrix_rho_threshold;
                warp_cell{i1,i2,i3} = nanmedian(cell2mat(warp_cell(i1,i2,(i3-half_medfilt_length-1)+find(correlated_id))),3);  
            end
        end
    end
end

% Clean warp_cell with a threshold for sudden matrix shifts
for i1 = 1:warp_blocksize
    for i2 = 1:warp_blocksize
        temp = find(sum(sum(abs(diff(cell2mat(warp_cell(i1,i2,:)),1,3)),1),2) > affinematrix_abssum_jump_threshold);
        if ~isempty(temp)
            for i3 = 1:length(temp)
                if strcmp(transform, 'translation')
                    warp_cell{i1,i2,temp(i3)} = NaN(2,1);
                    warp_cell{i1,i2,temp(i3)+1} = NaN(2,1);
                elseif strcmp(transform, 'homography')
                    warp_cell{i1,i2,temp(i3)} = NaN(3,3);
                    warp_cell{i1,i2,temp(i3)+1} = NaN(3,3);
                else
                    warp_cell{i1,i2,temp(i3)} = NaN(2,3);
                    warp_cell{i1,i2,temp(i3)+1} = NaN(2,3);
                end
            end
        end
    end
end
for i3 = 1:nz
    for i1 = 1:warp_blocksize
        for i2 = 1:warp_blocksize
            if isempty(warp_cell{i1,i2,i3})
                if strcmp(transform, 'translation')
                    warp_cell{i1,i2,i3} = NaN(2,1);
                elseif strcmp(transform, 'homography')
                    warp_cell{i1,i2,i3} = NaN(3,3);
                else
                    warp_cell{i1,i2,i3} = NaN(2,3);
                end
            end
        end
    end
end
for i1 = 1:warp_blocksize
    for i2 = 1:warp_blocksize
        temp = cell2mat(warp_cell(i1,i2,:));
        try
            for i3_1 = 1:2
                for i3_2 = 1:3
                    temp(i3_1,i3_2,:) = interp1(find(~isnan(temp(i3_1,i3_2,:))),permute(temp(i3_1,i3_2,~isnan(temp(i3_1,i3_2,:))),[3,1,2]),1:nz,'linear','extrap');
                end
            end
        catch
            temp = repmat(warp_init,[1,1,nz]);
        end
        for i3 = 1:nz
            warp_cell{i1,i2,i3} = temp(:,:,i3);
        end
    end
end

save(fn_affinetransmat, 'warp_cell', 'rho', 'qN_x', 'qN_y', 'nx', 'ny', 'nz', 'nonzero_row', 'nonzero_column',...
    'warp_template_tiffstack_num', 'warp_movave_tiffstack_num', 'warp_blocksize',...
    'warp_overlap_pix', 'n_split4warpinit', 'edge_remove_pix', 'affinematrix_abssum_threshold',...
    'affinematrix_rho_threshold', 'affinematrix_medfilt_tiffstack_num', 'downsample_frame_num');

%% Apply warp to tiff files
disp('Applying affine transformations...')
fns_summary = dir(fullfile(source_path, '*summary.mat'));
fns_tiff = dir(fullfile(source_path, '*.tif'));
if length(fns_summary) ~= length(fns_tiff)
   error('tiff number does not match summary file number');
end

if length(fns_tiff) ~= nz
   error('warp number does not match file number');
end

% apply the warp to each file and save
parfor i = 1:length(fns_tiff)
    applywarp_Npatches(fns_tiff(i).name, fns_summary(i).name, source_path, save_path, align_ch, n_ch, warp_cell(:,:,i), transform, edge_remove_pix, nonzero_row, nonzero_column, qN_x, qN_y, warp_overlap_pix, warp_blocksize);
end

%% Generate downsampled movie using summary file
disp('Generateing downsampled movies...')
fns_summary_warp = dir(fullfile(save_path, '*summary_warped.mat'));
if length(fns_summary_warp) ~= nz
   error('warp number does not match file number');
end

% load all summary files
% summary_temp = load(fullfile(save_path, fns_summary_warp(nz).name));
% downsampled_frame_num_perstack_last = size(summary_temp.downsampled, 3);
% summary_temp = load(fullfile(save_path, fns_summary_warp(nz-1).name));
% downsampled_frame_num_perstack = size(summary_temp.downsampled, 3);
% 
% downsampled = zeros(ny, nx, (nz-1)*downsampled_frame_num_perstack+downsampled_frame_num_perstack_last);
% downsampled_perstack = zeros(ny, nx, nz);
% for i = 1:nz-1
%     summary_temp = load(fullfile(save_path, fns_summary_warp(i).name));
%     downsampled(:, :, (i-1)*downsampled_frame_num_perstack+1:i*downsampled_frame_num_perstack) = int16(summary_temp.downsampled); 
%     downsampled_perstack(:, :, i) = int16(summary_temp.downsampled_perstack); 
% end
% summary_temp = load(fullfile(save_path, fns_summary_warp(nz).name));
% downsampled(:, :, (nz-1)*downsampled_frame_num_perstack+1:end) = int16(summary_temp.downsampled); 
% downsampled_perstack(:, :, nz) = int16(summary_temp.downsampled_perstack); 
% downsampled = int16(downsampled);
% downsampled_perstack = int16(downsampled_perstack);

downsampled_c = cell(1, 1, nz);
downsampled_perstack_c = cell(1, 1, nz);
parfor i = 1:nz
    summary_temp = load(fullfile(save_path, fns_summary_warp(i).name));
    downsampled_c{i} = int16(summary_temp.downsampled); 
    downsampled_perstack_c{i} = int16(summary_temp.downsampled_perstack); 
end
downsampled = cell2mat(downsampled_c);
downsampled_perstack = cell2mat(downsampled_perstack_c);

% save
write_tiff(fn_downsampled_perstack_save, downsampled_perstack, info_downsampled_perstack);
write_tiff(fn_downsampled2save, downsampled, info_downsampled_perstack);
write_tiff(fn_downsampled2save_max, max(downsampled, [], 3), info_downsampled_perstack);
write_tiff(fn_downsampled2save_mean, mean(downsampled, 3), info_downsampled_perstack);
%     write_tiff(fn_downsampled2save_max_norm, imnormalize2(max(downsampled, [], 3), affine_norm_radius));
%     write_tiff(fn_downsampled2save_mean_norm, imnormalize2(mean(downsampled, 3), affine_norm_radius));

disp('Warping complete')

