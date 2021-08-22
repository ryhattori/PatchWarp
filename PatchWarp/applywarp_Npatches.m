function applywarp_Npatches(fns_tiff_name, fns_summary_name, source_path, save_path, ch, n_ch, warp_cell, transform, edge_remove_pix, nonzero_row, nonzero_column, qN_x, qN_y, overlap_npix, block_size)
    % check if warped, summary file has been changed by adding a 'warp'
    % variable
    temp_summary = load(fullfile(source_path,fns_summary_name));
    if isfield(temp_summary,'warp')
        %disp(['skip warped file #',num2str(i)]); continue;
    end

    % if on local computer, copy image to temporary file, otherwise load
%     if ~ispc
        img_filename = fullfile(source_path, fns_tiff_name);
%     else
% %         disp('Copying tiff to local drive...')
%         img_filename = fullfile(source_path, fns_tiff_name);
% %         disp(img_filename);
%         tic
%         copyfile(img_filename,local_filename);
%         toc
% %         disp('Done.')
%         img_filename = local_filename;
% %         disp(img_filename);
%     end
    
    [stack,info, ~] = read_tiff(img_filename, ch, n_ch);

    stack = stack(:,1+edge_remove_pix:end-edge_remove_pix,:); % remove n pixels from left and right edges.
    if size(stack,1) == length(nonzero_row) && size(stack,2) == length(nonzero_column)
        stack = stack(nonzero_row,nonzero_column,:); % remove zero-zone
    end
    stack_warp = nan(size(stack,1),size(stack,2),size(stack,3));
    temp_summary.warp_cell = warp_cell;
    %apply warping to each frame
    for i3 = 1:size(stack,3)
        out_temp = cell(block_size,block_size);
        for i1 = 1:block_size
            for i2 = 1:block_size
                out_temp{i1,i2} = spatial_interp_patchwarp(double(stack(qN_y{i1,i2},qN_x{i1,i2},i3)), warp_cell{i1,i2}, transform, 1:length(qN_x{i1,i2}), 1:length(qN_y{i1,i2}));
                out_temp{i1,i2}(out_temp{i1,i2}==0) = NaN;
            end
        end
        
        out_temp = cell2mat(out_temp);
        for i1 = 2:block_size
            out_temp(qN_y{i1-1,1}(end)-(2*overlap_npix+1)+1:qN_y{i1-1,1}(end),:) = ...
                nanmean(cat(3,out_temp(qN_y{i1-1,1}(end)-(2*overlap_npix+1)+1:qN_y{i1-1,1}(end),:)...
                , out_temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*overlap_npix+1),:)),3);
            out_temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*overlap_npix+1),:) = [];
        end
        for i2 = 2:block_size
            out_temp(:,qN_x{1,i2-1}(end)-(2*overlap_npix+1)+1:qN_x{1,i2-1}(end)) = ...
                nanmean(cat(3,out_temp(:,qN_x{1,i2-1}(end)-(2*overlap_npix+1)+1:qN_x{1,i2-1}(end))...
                , out_temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*overlap_npix+1))),3);
            out_temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*overlap_npix+1)) = [];
        end
%         out_temp(isnan(out_temp))=0;
        nan_elements = find(isnan(out_temp));
        out_temp(nan_elements) = (prctile(out_temp(~isnan(out_temp)),30)-prctile(out_temp(~isnan(out_temp)),1)).*rand(length(nan_elements),1) + prctile(out_temp(~isnan(out_temp)),1);

        stack_warp(:,:,i3) = int16(out_temp);
    end
    % save warped file, replace local
    % not local PC: first save a temp file in folder, then replace the file
    % before warp
    % replace summary file
    if ~ispc
%         disp(['save temp file:', fns_tiff_name(1:end-4),'_warped.tif']);
        write_tiff(fullfile(save_path, [fns_tiff_name(1:end-4),'_warped.tif']),...
            int16(stack_warp),info);
        %disp('replace original tiff')
        %[~,msg] = copyfile(fullfile(source_path, [fns_tiff_name(1:end-4),'_warped.tif']), ...
        %    fullfile(source_path, 'motioncorrected_tiff',fns_tiff_name));
        %disp(msg);
        %delete
        %delete(fullfile(source_path, [fns_tiff_name(1:end-4),'_warped.tif']));
        %replace summary file
%         disp('replace summry file')
        save(fullfile(save_path,[fns_summary_name(1:end-4) '_warped.mat']), '-struct','temp_summary');
        clear temp_summary
    else
        %replace local temp file and copy to server
        write_tiff(img_filename,...
            int16(stack_warp),info);
%         disp('replace original tiff')
        [~,msg] = copyfile(img_filename, ...
            fullfile(save_path,[fns_tiff_name(1:end-4) '_warped.tif']));
%         disp(msg);
%         disp('replace summry file')
        save(fullfile(save_path,[fns_summary_name(1:end-4) '_warped.mat']), '-struct','temp_summary');
        clear temp_summary
    end
end
