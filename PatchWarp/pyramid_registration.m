function done = pyramid_registration(fn, target, save_path, align_ch, save_ch, n_downsampled, n_downsampled_perstack, n_ch, rigid_norm_method, rigid_norm_radius, rigid_template_center_frac)
    % save_path should be made beforehand, and all the arguments should be given.
    % if target is empty, do nothing but return whether it is done.
    % 
    % Original code was written by Aki Mitani and available at https://github.com/amitani/matlab_motion_correct
    % revised by Ryoma Hattori for PatchWarp pipeline
    %%
    
    [~,fn_root]=fileparts(fn);
    fn_corrected_tif = fullfile(save_path,[fn_root '_corrected.tif']);
    fn_summary_mat = fullfile(save_path,[fn_root '_summary.mat']);
    
    if(isempty(target))
        done = java.io.File(fn_corrected_tif).exists() ...
               & java.io.File(fn_summary_mat).exists();
        return;
    end
    
    L = Logger();
    
    if(ischar(target))
        L.newline('Reading target image. %s',target);
    else
        L.newline('Reading target image.');
    end
    target=parse_image_input(target,align_ch);
    
%     if(isempty(n_ch))
%         L.newline('Done. Checking the number of channels. %s', fn);
%         si4 = read_SI4_header(fn);
%         n_ch = numel(si4.channelsSave);
%     end
    
    L.newline('Done. Reading source image. %s', fn);
    if length(save_ch) > 1
        [image_stack_align, ~] = read_tiff(fn,align_ch,n_ch);
        [image_stack_save, info] = read_tiff(fn,1,1);
    else
        if(align_ch == save_ch)
            [image_stack_save, info] = read_tiff(fn,save_ch,n_ch);
            image_stack_align = image_stack_save;
        else
            [image_stack, info] = read_tiff(fn,[align_ch save_ch],n_ch);
            image_stack_align = image_stack(:,:,:,1);
            image_stack_save = image_stack(:,:,:,2);
        end
    end
    
    if strcmp(rigid_norm_method, 'rank')
        target = rank_transform(target);
        image_stack_align = rank_transform(image_stack_align);
    elseif strcmp(rigid_norm_method, 'local')
        target = imnormalize2(double(target), rigid_norm_radius);
        image_stack_align = imnormalize2(double(image_stack_align), rigid_norm_radius);
    else
        disp(['Warning: ' rigid_norm_method 'is not a valid normalization method!!'])
    end
    
    L.newline('Done. Motion correcting.');
    if length(save_ch) > 1
        ir = BilinearPyramidImageRegistrator(target, rigid_template_center_frac, 3);
        t = zeros(size(image_stack_save,3),2);
        for j = 1:size(image_stack_align,3)
            t(2*j-1,:) = ir.register(double(image_stack_align(:,:,j)));
            t(2*j,:) = t(2*j-1,:);
        end
    else
        ir = BilinearPyramidImageRegistrator(target, rigid_template_center_frac, 3);
        t = zeros(size(image_stack_align,3),2);
        for j = 1:size(image_stack_align,3)
            t(j,:)=ir.register(double(image_stack_align(:,:,j)));
        end
    end
    
    L.newline('Done. Shifting signal ch.');
    corrected = zeros(size(image_stack_save),'single');
    for j = 1:size(image_stack_save,3)
        corrected(:,:,j)=BilinearPyramidImageRegistrator.shift(...
            image_stack_save(:,:,j),t(j,:));
    end

    
%     % Added by RH to change scanimage.SI4.channelsSave from original to save_ch
%     if ~isequal(str2double(regexp(info(1).ImageDescription,'channelsSave = (\d*)','tokens','once')),save_ch)
%         for i = 1:numel(info)
%             [~,endIndex1] = regexp(info(i).ImageDescription,'channelsSave = ');
%             [~,endIndex2] = regexp(info(i).ImageDescription,['channelsSave = (\d*)' newline]);
%             info(i).ImageDescription(endIndex1+1:endIndex2-1) = num2str(save_ch);
%         end
%     end
%     %
    
    L.newline('Done. Saving corrected files.');
    write_tiff(fn_corrected_tif,int16(corrected),info);
    
    L.newline('Done. Making downsampled movie.');
    if(n_downsampled>0)
        downsampled = cast(downsample_mean(corrected,n_downsampled,3),class(image_stack_save));
    else
        downsampled = [];
    end
        
    if(n_downsampled_perstack>0)
        L.newline('Done. Shifting align ch.');
        if length(save_ch) > 1
            corrected_align = zeros(size(image_stack_align),'single');
            for j = 1:size(image_stack_align,3)
                corrected_align(:,:,j)=BilinearPyramidImageRegistrator.shift(...
                    image_stack_align(:,:,j),t(2*j-1,:));
            end
        else
            corrected_align = zeros(size(image_stack_align),'single');
            for j = 1:size(image_stack_align,3)
                corrected_align(:,:,j)=BilinearPyramidImageRegistrator.shift(...
                    image_stack_align(:,:,j),t(j,:));
            end
        end
        downsampled_perstack = cast(downsample_mean(corrected_align, n_downsampled_perstack,3), class(image_stack_align));
    else
        downsampled_perstack = [];
    end
    
    
    L.newline('Done. Saving data.');
    
%     [frame_tag, si4, info_first] = parse_SI4_info(info);
    
    method = 'interpolate';
    
    if strcmp(rigid_norm_method, 'rank')
        method = [method ' rank']; 
    elseif strcmp(rigid_norm_method, 'local')
        method = [method ' localnorm']; 
    end

    save(fn_summary_mat, '-v6',   'downsampled',...
                                  'downsampled_perstack',...
                                  'info',...
                                  't',...
                                  'target',...
                                  'method', ...
                                  'align_ch',...
                                  'save_ch',...
                                  'n_downsampled',...
                                  'n_downsampled_perstack');


    L.newline('Done.');
    
    done=true;
end

