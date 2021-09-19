function [input_images_warp1, input_images_warp2] = patchwarp_across_sessions_apply(input_images, patchwarp_results, second_affine)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PatchWarp application to registration between different imaging sessions
% -------------------
% 
% Released by Ryoma Hattori
% Email: rhattori0204@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Adjust the input image size
if size(input_images, 2) < patchwarp_results.xmax
    input_images = cat(2, input_images, zeros(size(input_images, 1), patchwarp_results.xmax - size(input_images, 2), size(input_images, 3)));
%     input_images = cat(2, input_images, (prctile(input_images(:), 30) - prctile(input_images(:), 1)).*rand(size(input_images, 1), patchwarp_results.xmax - size(input_images, 2), size(input_images, 3)) + prctile(input_images(:), 1));
end
if size(input_images, 1) < patchwarp_results.ymax
    input_images = cat(1, input_images, zeros(patchwarp_results.ymax - size(input_images, 1), size(input_images, 2), size(input_images, 3)));
end
    
%% Apply the 1st affine transformation to the whole FOV
input_images_warp1 = zeros(size(input_images, 1), size(input_images, 2), size(input_images, 3));
[~, rho1_max_id] = max(patchwarp_results.rho1);
for i = 1:size(input_images, 3)
    input_images_warp1(:, :, i) = spatial_interp_patchwarp(input_images(:, :, i), patchwarp_results.warp1_cell{rho1_max_id}, patchwarp_results.transform1, 1:patchwarp_results.xmax, 1:patchwarp_results.ymax);
end

%% Apply the 2nd affine transformations to the subfields if second_affine==1
if second_affine==1
    input_images_warp2 = zeros(size(input_images, 1), size(input_images, 2), size(input_images, 3));
    for i = 1:size(input_images, 3)
        input_images_warp2_temp = cell(patchwarp_results.warp_blocksize, patchwarp_results.warp_blocksize);
        for i1 = 1:patchwarp_results.warp_blocksize
            for i2 = 1:patchwarp_results.warp_blocksize
                [~, rho2_max_id] = max(patchwarp_results.rho2(i1, i2, :));
                input_images_warp2_temp{i1, i2} = spatial_interp_patchwarp(double(input_images_warp1(patchwarp_results.qN_y{i1, i2}, patchwarp_results.qN_x{i1, i2}, i)), patchwarp_results.warp2_cell{i1, i2, rho2_max_id}, patchwarp_results.transform2, 1:length(patchwarp_results.qN_x{i1, i2}), 1:length(patchwarp_results.qN_y{i1, i2}));
                input_images_warp2_temp{i1, i2}(input_images_warp2_temp{i1, i2}==0) = NaN;
            end
        end
        input_images_warp2_temp = cell2mat(input_images_warp2_temp);
        for i1 = 2:patchwarp_results.warp_blocksize
            input_images_warp2_temp(patchwarp_results.qN_y{i1-1, 1}(end)-(2*patchwarp_results.warp_overlap_pix +1)+1:patchwarp_results.qN_y{i1-1, 1}(end), :) = ...
                nanmean(cat(3,input_images_warp2_temp(patchwarp_results.qN_y{i1-1, 1}(end)-(2*patchwarp_results.warp_overlap_pix +1)+1:patchwarp_results.qN_y{i1-1, 1}(end), :)...
                , input_images_warp2_temp(patchwarp_results.qN_y{i1-1, 1}(end)+1:patchwarp_results.qN_y{i1-1, 1}(end)+(2*patchwarp_results.warp_overlap_pix +1), :)),3);
            input_images_warp2_temp(patchwarp_results.qN_y{i1-1, 1}(end)+1:patchwarp_results.qN_y{i1-1, 1}(end)+(2*patchwarp_results.warp_overlap_pix +1), :) = [];
        end
        for i2 = 2:patchwarp_results.warp_blocksize
            input_images_warp2_temp(:, patchwarp_results.qN_x{1, i2-1}(end)-(2*patchwarp_results.warp_overlap_pix +1)+1:patchwarp_results.qN_x{1, i2-1}(end)) = ...
                nanmean(cat(3, input_images_warp2_temp(:, patchwarp_results.qN_x{1, i2-1}(end)-(2*patchwarp_results.warp_overlap_pix +1)+1:patchwarp_results.qN_x{1, i2-1}(end))...
                , input_images_warp2_temp(:, patchwarp_results.qN_x{1, i2-1}(end)+1:patchwarp_results.qN_x{1, i2-1}(end)+(2*patchwarp_results.warp_overlap_pix +1))), 3);
            input_images_warp2_temp(:, patchwarp_results.qN_x{1, i2-1}(end)+1:patchwarp_results.qN_x{1, i2-1}(end)+(2*patchwarp_results.warp_overlap_pix +1)) = [];
        end
        nan_elements = find(isnan(input_images_warp2_temp));
        input_images_warp2_temp(nan_elements) = zeros(length(nan_elements), 1);
        input_images_warp2(:, :, i) = input_images_warp2_temp;
    end
else
    input_images_warp2 = [];
end

end