function applywarp_Npatches_singletiffstack(fns_tiff_name, fns_summary_name, source_path, save_path, ch, n_ch, warp_cell, transform, edge_remove_pix, nonzero_row, nonzero_column, qN_x, qN_y, overlap_npix, block_size)
    temp_summary = load(fullfile(source_path, fns_summary_name));
    img_filename = fullfile(source_path, fns_tiff_name);
    
    [stack, info, ~] = read_tiff(img_filename, ch, n_ch);

    stack = stack(:, 1+edge_remove_pix:end-edge_remove_pix, :); % remove n pixels from left and right edges.
    if size(stack, 1) == length(nonzero_row) && size(stack, 2) == length(nonzero_column)
        stack = stack(nonzero_row, nonzero_column, :); % remove zero-zone
    end
    
    stack_warp = nan(length(nonzero_row), length(nonzero_column), size(stack, 3));

    temp_summary.warp_cell = warp_cell;
    %apply warping to each frame
    parfor i3 = 1:size(stack, 3)
        out_temp = cell(block_size, block_size);
        for i1 = 1:block_size
            for i2 = 1:block_size
                out_temp{i1, i2} = spatial_interp_patchwarp(double(stack(qN_y{i1,i2}, qN_x{i1,i2}, i3)), warp_cell{i1, i2, ceil(i3/temp_summary.n_downsampled)}, transform, 1:length(qN_x{i1, i2}), 1:length(qN_y{i1, i2}));
                out_temp{i1, i2}(out_temp{i1, i2}==0) = NaN;
            end
        end
        
        out_temp = cell2mat(out_temp);
        for i1 = 2:block_size
            out_temp(qN_y{i1-1, 1}(end)-(2*overlap_npix+1)+1:qN_y{i1-1, 1}(end), :) = ...
                nanmean(cat(3,out_temp(qN_y{i1-1, 1}(end)-(2*overlap_npix+1)+1:qN_y{i1-1, 1}(end), :)...
                , out_temp(qN_y{i1-1, 1}(end)+1:qN_y{i1-1, 1}(end)+(2*overlap_npix+1), :)),3);
            out_temp(qN_y{i1-1, 1}(end)+1:qN_y{i1-1, 1}(end)+(2*overlap_npix+1), :) = [];
        end
        for i2 = 2:block_size
            out_temp(:, qN_x{1, i2-1}(end)-(2*overlap_npix+1)+1:qN_x{1, i2-1}(end)) = ...
                nanmean(cat(3, out_temp(:, qN_x{1, i2-1}(end)-(2*overlap_npix+1)+1:qN_x{1, i2-1}(end))...
                , out_temp(:, qN_x{1, i2-1}(end)+1:qN_x{1, i2-1}(end)+(2*overlap_npix+1))), 3);
            out_temp(:, qN_x{1, i2-1}(end)+1:qN_x{1, i2-1}(end)+(2*overlap_npix+1)) = [];
        end
        
        % Fill nan with either 0 or background level
%         out_temp(isnan(out_temp))=0;
        out_temp2 = nan(length(nonzero_row), length(nonzero_column));
        out_temp2(nonzero_row, nonzero_column) = out_temp;
        out_temp = out_temp2;
%         nan_elements = find(isnan(out_temp));
%         out_temp(nan_elements) = (prctile(out_temp(~isnan(out_temp)), 30) - prctile(out_temp(~isnan(out_temp)), 1)).*rand(length(nan_elements), 1) + prctile(out_temp(~isnan(out_temp)), 1);

        stack_warp(:, :, i3) = int16(out_temp);
    end
    
    % Downsample warped images
    temp_summary.downsampled_perstack = downsample_mean(stack_warp, temp_summary.n_downsampled_perstack, 3);
    temp_summary.downsampled = downsample_mean(stack_warp, temp_summary.n_downsampled, 3);
    
    % Save
    write_tiff(fullfile(save_path, [fns_tiff_name(1:end-4), '_warped.tif']), int16(stack_warp), info);
    save(fullfile(save_path,[fns_summary_name(1:end-4) '_warped.mat']), '-struct', 'temp_summary');
    clear temp_summary
end
