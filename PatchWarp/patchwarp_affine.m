function patchwarp_affine(source_path, save_path, n_ch, align_ch, warp_template_tiffstack_num, warp_movave_tiffstack_num, warp_blocksize, warp_overlap_pix, n_split4warpinit, edge_remove_pix, affinematrix_abssum_threshold, affinematrix_rho_threshold, affinematrix_medfilt_tiffstack_num, downsample_frame_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PatchWarp
% pacman_warp
% -------------------
% Piece-wise affine transformation of rigid-motion corrected images
% 
% Released by Ryoma Hattori
% Email: rhattori0204@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
   %make temp folder and files
   local_dir = 'C:\temp_files';
   local_filename = [local_dir filesep 'temp_img_' datestr(now,'yymmddHHMMSSFFF') '.tif'];
   if ~isdir(local_dir)
       mkdir(local_dir)
   end
   % delete local temp file if it exists
   if exist(local_filename,'file')
       delete(local_filename)
   end
end

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
fn_downsampled = fullfile(source_path,'downsampled',fn_downsampled_name);
fn_downsampled2save = fullfile(save_path,'downsampled',[fn_downsampled_name(1:end-4),'_warped.tif']);
fn_downsampled2save_max = fullfile(save_path,'downsampled',[fn_downsampled_name(1:end-4),'_warped_max.tif']);
fn_downsampled2save_mean = fullfile(save_path,'downsampled',[fn_downsampled_name(1:end-4),'_warped_mean.tif']);
% fn_downsampled2save_max_norm = fullfile(save_path,'downsampled',[fn_downsampled_name(1:end-4),'_warped_max_normalized.tif']);
% fn_downsampled2save_mean_norm = fullfile(save_path,'downsampled',[fn_downsampled_name(1:end-4),'_warped_mean_normalized.tif']);

%default parameters for warping
% n_ch = 1;
% ch = 1;
% margin_ratio = 0.125;
n_depth = 1;
n_step = 50;
transform = 'affine';
normalize_r1 = 1;
%noramlize_r2 = 10; % default 32. Decrease this number if image quality is bad (dark). For example, 10.
 noramlize_r2 = 32; % Use when you use non-shrunk images for estimating
% affine transformation?
noramlize_offset = 50;

%load tif files
disp('load downsampled_perstack file')

if ispc %copy to local
    disp('Copying downsampled_perstack to local drive...')
    img_filename = fn_downsampled_perstack;
    disp(img_filename);
    copyfile(img_filename,local_filename);
    disp('Done.')
    img_filename = local_filename;
    disp(img_filename);
    [stack_downsampled_perstack, info_downsampled_perstack] = read_tiff(img_filename, align_ch, n_ch);
else
    [stack_downsampled_perstack, info_downsampled_perstack] = read_tiff(fn_downsampled_perstack, align_ch, n_ch);
end

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
qN_x = cell(warp_blocksize,warp_blocksize);
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

% first get the warp transform, then median filter it to use
disp('warping downsampled_perstack image')

%get warp transform
warp_cell = cell(warp_blocksize, warp_blocksize, nz);

if mod(n_split4warpinit, 2) == 1
    error('Session split must be even number for warping')
end

template_im = mean(stack_downsampled_perstack(:, :, round(nz/2)-round((warp_template_tiffstack_num - 1)/2):round(nz/2)- round((warp_template_tiffstack_num - 1)/2) + warp_template_tiffstack_num - 1), 3); % default (:,:,1)
template_im_normalized = imnormalize(template_im, normalize_r1, noramlize_r2, noramlize_offset);

smoothed_downsampled_perstack = NaN(ny, nx, nz);
for i = 1:nz
    if i <= round((warp_movave_tiffstack_num - 1) / 2)
        smoothed_downsampled_perstack(:,:,i) = int16(imnormalize(mean(stack_downsampled_perstack(:, :, 2:warp_movave_tiffstack_num), 3), normalize_r1, noramlize_r2, noramlize_offset));
    elseif i > nz - round((warp_movave_tiffstack_num - 1) / 2)
        smoothed_downsampled_perstack(:,:,i) = int16(imnormalize(mean(stack_downsampled_perstack(:, :, nz-warp_movave_tiffstack_num:nz-1), 3), normalize_r1, noramlize_r2, noramlize_offset));
    else
        smoothed_downsampled_perstack(:,:,i) = int16(imnormalize(mean(stack_downsampled_perstack(:, :, i-round((warp_movave_tiffstack_num - 1) / 2):i-round((warp_movave_tiffstack_num - 1) / 2) + warp_movave_tiffstack_num - 1),3), normalize_r1, noramlize_r2, noramlize_offset));
    end
end

affine_init = [1 0 0; 0 1 0];
warp4template = cell(warp_blocksize, warp_blocksize, n_split4warpinit);
for i1 = 1:warp_blocksize
    for i2 = 1:warp_blocksize
        warp4template{i1,i2,n_split4warpinit/2} = affine_init;
        warp4template{i1,i2,1+n_split4warpinit/2} = affine_init;
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

n_warp2use4template = 7;
for nth_warp = [n_split4warpinit/2:-1:1,n_split4warpinit/2+1:n_split4warpinit]
    for i3 = warp_range(nth_warp,1):warp_range(nth_warp,2)
        parfor i1 = 1:warp_blocksize
            for i2 = 1:warp_blocksize
                  [~, warp_cell{i1,i2,i3}, ~, rho(i1,i2,i3)]= ...
                    ecc_patchwarp(smoothed_downsampled_perstack(qN_y{i1,i2}, qN_x{i1,i2},i3), template_im_normalized(qN_y{i1,i2}, qN_x{i1,i2}), n_depth, n_step, transform, warp4template{i1, i2, nth_warp}(1:2,1:3));
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
                    warp_cell{i1,i2,i3} = NaN(2,3);
                end
                if warp_cell{i1,i2,i3}(1,1)<0.6 || warp_cell{i1,i2,i3}(2,2)<0.6 || warp_cell{i1,i2,i3}(1,1)>1.4 || warp_cell{i1,i2,i3}(2,2)>1.4 || rho(i1,i2,i3) <= affinematrix_rho_threshold
                    warp_cell{i1,i2,i3} = NaN(2,3);
                end
                if size(warp_cell{i1,i2,i3},1) == 3
                    warp_cell{i1,i2,i3}(3,:) = [];
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

SuddenJumpThreshold = 10;
for i1 = 1:warp_blocksize
    for i2 = 1:warp_blocksize
        temp = find(sum(sum(abs(diff(cell2mat(warp_cell(i1,i2,:)),1,3)),1),2) > SuddenJumpThreshold);
        if ~isempty(temp)
            for i3 = 1:length(temp)
                warp_cell{i1,i2,temp(i3)} = NaN(2,3);
                warp_cell{i1,i2,temp(i3)+1} = NaN(2,3);
            end
        end
    end
end

for i3 = 1:nz
    for i1 = 1:warp_blocksize
        for i2 = 1:warp_blocksize
            if isempty(warp_cell{i1,i2,i3})
                warp_cell{i1,i2,i3} = NaN(2,3);
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
            temp = repmat(affine_init,[1,1,nz]);
        end
        for i3 = 1:nz
            warp_cell{i1,i2,i3} = temp(:,:,i3);
        end
    end
end

% apply warp
out_downsampled_perstack = cell(warp_blocksize, warp_blocksize, nz);
for i3 = 1:nz
    for i1 = 1:warp_blocksize
        for i2 = 1:warp_blocksize
            out_downsampled_perstack{i1,i2,i3} = spatial_interp_patchwarp(double(stack_downsampled_perstack(qN_y{i1,i2}, qN_x{i1,i2},i3)), warp_cell{i1,i2,i3}, transform, 1:length(qN_x{i1,i2}), 1:length(qN_y{i1,i2}));
            out_downsampled_perstack{i1,i2,i3}(out_downsampled_perstack{i1,i2,i3}==0) = NaN;
        end
    end
end

stack_downsampled_perstack_warp = NaN(ny, nx, nz);
for i3 = 1:nz
    temp = cell2mat(out_downsampled_perstack(:,:,i3));
    for i1 = 2:warp_blocksize
        temp(qN_y{i1-1,1}(end)-(2*warp_overlap_pix+1)+1:qN_y{i1-1,1}(end),:) = ...
            nanmean(cat(3,temp(qN_y{i1-1,1}(end)-(2*warp_overlap_pix+1)+1:qN_y{i1-1,1}(end),:)...
            , temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*warp_overlap_pix+1),:)),3);
        temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*warp_overlap_pix+1),:) = [];
    end
    for i2 = 2:warp_blocksize
        temp(:,qN_x{1,i2-1}(end)-(2*warp_overlap_pix+1)+1:qN_x{1,i2-1}(end)) = ...
            nanmean(cat(3,temp(:,qN_x{1,i2-1}(end)-(2*warp_overlap_pix+1)+1:qN_x{1,i2-1}(end))...
            , temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*warp_overlap_pix+1))),3);
        temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*warp_overlap_pix+1)) = [];
    end
    
    temp(isnan(temp))=0;
    stack_downsampled_perstack_warp(:,:,i3) = int16(temp); 
end

% save downsampled_perstack file and copy to server (if local)
disp('save wapred downsampled_perstack file and affine transformation matrices')

write_tiff(fn_downsampled_perstack_save, int16(stack_downsampled_perstack_warp), info_downsampled_perstack);
save(fn_affinetransmat, 'warp_cell', 'rho');


%% generate downsampled movie using summary file; summary file is fast
% to load no need for local copy. 
disp('Warp and generate downsampled movie')
fns_summary = dir(fullfile(source_path, '*summary.mat'));
if length(fns_summary) ~= nz
   error('warp number does not match file number');
end
if any(exist(fn_downsampled2save, 'file')) % check if done
    disp('Warped downsampled movie has been generated, skip')
else
    % load all summary files
    downsampled_c = cell(1,1,nz);
    for i=1:nz
        summary_temp=load(fullfile(source_path, fns_summary(i).name));
        downsampled_c{i} = summary_temp.downsampled(:,1+edge_remove_pix:end-edge_remove_pix, :); % removed n pix from left and right. 
        downsampled_c{i} = downsampled_c{i}(nonzero_row, nonzero_column,:); % remove zero-zone
    end

    % warp from the 2nd to the last
    disp('Apply warp') % NOTE: summary files is not modified by warp

    parfor i3 = 1:nz
         %apply warping to each frame
        for i_fr = 1:size(downsampled_c{i3}, 3)
            out_temp = cell(warp_blocksize, warp_blocksize);
            for i1 = 1:warp_blocksize
                for i2 = 1:warp_blocksize
                    out_temp{i1, i2} = spatial_interp_patchwarp(double(downsampled_c{i3}(qN_y{i1, i2}, qN_x{i1, i2}, i_fr)), warp_cell{i1, i2, i3}, transform, 1:length(qN_x{i1, i2}), 1:length(qN_y{i1, i2}));
                    out_temp{i1, i2}(out_temp{i1, i2}==0) = NaN;
                end
            end

            out_temp = cell2mat(out_temp);
            for i1 = 2:warp_blocksize
                out_temp(qN_y{i1-1,1}(end)-(2*warp_overlap_pix+1)+1:qN_y{i1-1,1}(end),:) = ...
                    nanmean(cat(3,out_temp(qN_y{i1-1,1}(end)-(2*warp_overlap_pix+1)+1:qN_y{i1-1,1}(end),:)...
                    , out_temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*warp_overlap_pix+1),:)),3);
                out_temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*warp_overlap_pix+1),:) = [];
            end
            for i2 = 2:warp_blocksize
                out_temp(:,qN_x{1,i2-1}(end)-(2*warp_overlap_pix+1)+1:qN_x{1,i2-1}(end)) = ...
                    nanmean(cat(3,out_temp(:,qN_x{1,i2-1}(end)-(2*warp_overlap_pix+1)+1:qN_x{1,i2-1}(end))...
                    , out_temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*warp_overlap_pix+1))),3);
                out_temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*warp_overlap_pix+1)) = [];
            end

            out_temp(isnan(out_temp)) = 0;

            downsampled_c{i3}(:,:,i_fr) = int16(out_temp); 
        end
    end
    downsampled = cell2mat(downsampled_c);

    %save
    disp('save warped movie')
    write_tiff(fn_downsampled2save, downsampled);
    write_tiff(fn_downsampled2save_max, max(downsampled, [], 3));
    write_tiff(fn_downsampled2save_mean, mean(downsampled, 3));
%     write_tiff(fn_downsampled2save_max_norm, imnormalize(max(downsampled, [], 3), normalize_r1, noramlize_r2, noramlize_offset));
%     write_tiff(fn_downsampled2save_mean_norm, imnormalize(mean(downsampled, 3), normalize_r1, noramlize_r2, noramlize_offset));
end

%% apply warp to tiff files
fns_tiff = dir(fullfile(source_path, '*.tif'));
if length(fns_summary) ~= length(fns_tiff)
   error('tiff number does not match summary file number');
else

end
if length(fns_tiff) ~= nz
   error('warp number does not match file number');
end

%apply the warp to each file and save
parfor i = 1:length(fns_tiff)
    applywarp_Npatches(fns_tiff(i).name, fns_summary(i).name, source_path, save_path, align_ch, n_ch, warp_cell(:,:,i), transform, edge_remove_pix, nonzero_row, nonzero_column, qN_x, qN_y, warp_overlap_pix, warp_blocksize);
end

%% save warped downsampled movie to local folder if needed
% if nargin > 1
%     [~,local_downsampled_fn2save] = fileparts(fn_downsampled2save);
%     disp('save warped downsampled movie to local')
%     tic
%     copyfile(fn_downsampled2save,fullfile(save_path, [local_downsampled_fn2save,'.tif']));
%     toc
% end

disp('done warping')

