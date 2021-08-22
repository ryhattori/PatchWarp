function aligned = process_zstack_RH(fileNames)
%     fn = 'Z:\Data\ImagingRig5\180124\AM4170\AM4170_180124_002_001.tif';
%     fn = 'Z:\Data\ImagingRig3\171230\NH049\NH049_1230_00002_00001.tif';
%     fileNames = {};
    if ispc
        use_mex = 1;
    else
        use_mex = 0;
    end
use_mex=0
    if(iscell(fileNames))
        fn = fileNames{1};
    else
        fn = fileNames;
        fileNames = {};
    end
    %%
    info = imfinfo(fn);
    if(info(1).ImageDescription(1)=='f')
        %ScanImage 5
        SI = assignments2StructOrObj(info(1).Software);

        loggingFramesPerFile = SI.hScan2D.logFramesPerFile;
        channelsSave = SI.hChannels.channelSave;
        acqNumFrames = SI.hStackManager.framesPerSlice;
        stackNumSlices = SI.hStackManager.numSlices;
        fastZEnable = SI.hFastZ.enable;
        if(fastZEnable)
            acqNumFrames = SI.hFastZ.numVolumes;
            stackNumSlices  = SI.hFastZ.numFramesPerVolume;
        end % SI5 can do slow zstack with fastz, but this script is not compatible.
        clear SI
    else
        %ScanImage 4
        tmp = assignments2StructOrObj(info(1).ImageDescription);

        loggingFramesPerFile = tmp.SI4.loggingFramesPerFile;
        channelsSave = tmp.SI4.channelsSave;
        acqNumFrames = tmp.SI4.acqNumFrames;
        stackNumSlices = tmp.SI4.stackNumSlices;
        fastZEnable = tmp.SI4.fastZEnable;
        if(fastZEnable)
            acqNumFrames = tmp.SI4.fastZNumVolumes;
        end
        clear tmp;
    end
    numChannels = numel(channelsSave);
    %%
    totalFrames = acqNumFrames * stackNumSlices;
    numFiles = ceil(totalFrames / loggingFramesPerFile);

    if(isempty(fileNames))
        if(numFiles == 1)
            fileNames = {fn};
        else
            fileNames = cell(numFiles,1);
            for i=1:numFiles
                fileNames{i} =  [fn(1:end-7) sprintf('%03d',i) '.tif']; % assumes tif extension
            end
        end
    end
    disp(fileNames);
    %%
    im = cell(numFiles,1);
    

        for i=1:numFiles
            if exist(fileNames{i})
                disp(['Loading ' fileNames{i}]);
                im{i} = read_tiff(fileNames{i});
    %           im{i}(im{i}<80) = 35;
            else
            end
        end

%     dims = cell2mat(cellfun(@size, im(~cellfun(@isempty,im)), 'uni', false));
%     ydim = dims(1);
%     xdim = dims(2);
%     totaldims = sum(dims,3);
    im = cat(3,im{~cellfun(@isempty,im)});
    numFrames = size(im,3)/numChannels;
    im = reshape(im,size(im,1),size(im,2),numChannels,numFrames);
    %%
    if(fastZEnable)
        if(mod(numFrames,stackNumSlices))
            disp([numFrames acqNumFrames stackNumSlices]);
            warning('FastZ: # frames cannot be divided by stackNumSlices');
        end
        z = mod(0:numFrames-1, stackNumSlices)+1;
    else
        if(numFrames ~= acqNumFrames * stackNumSlices)
            disp([numFrames acqNumFrames stackNumSlices]);
            warning('NonFastZ: # frames doesn''t match');
        end
        z = floor((0:numFrames-1)/acqNumFrames)+1;
    end
    %%
    disp('Correcting each slice');
    stackNumSlices = max(z);
    correctedZStack = zeros(size(im,1),size(im,2),numChannels, stackNumSlices);
    ch_align = size(im,3);
    for i=1:stackNumSlices
        fprintf('%d/%d\n',i,stackNumSlices);
        im_slice = im(:,:,:,z==i);
        [correctedZStack(:,:,ch_align,i), t] = make_stable_average_RH(squeeze(im_slice(:,:,ch_align,:)), ch_align, use_mex);
        for ch=1:numChannels
            if(ch == ch_align)
                continue;
            end
            im_ch = zeros(size(im_slice,1),size(im_slice,2),size(im_slice,4));
            for j=1:size(im_slice,4)
                im_ch(:,:,j)=BilinearImageRegistrator.shift(im_slice(:,:,j),t(j,:));
            end
            correctedZStack(:,:,ch,i) = mean(im_ch,3);
        end
    end

    % StackViewer(correctedZStack)
    %%
    disp('Aligning slices');
    center_z_index = floor((stackNumSlices+1)/2);
    cumulative_t = cell(stackNumSlices,1);
    aligned = zeros(size(correctedZStack));
    aligned(:,:,:,center_z_index ) = correctedZStack(:,:,:,center_z_index);
    for direction = [-1 1]
        if(direction>0)
            max_dz=stackNumSlices-center_z_index;
        else
            max_dz=center_z_index-1;
        end
        cumulative_t{center_z_index}=[0 0];
        for dz = 1:max_dz
            if use_mex
                tz = cvMotionCorrect...
                    (correctedZStack(:,:,ch_align,center_z_index+direction*dz),...
                    correctedZStack(:,:,ch_align,center_z_index+direction*(dz-1)));
            else
                source = correctedZStack(:,:,ch_align,center_z_index+direction*dz);
                target = correctedZStack(:,:,ch_align,center_z_index+direction*(dz-1));
                ir = BilinearPyramidImageRegistrator(target,0.75,3);
                tz = zeros(size(source,3),2);
                for j = 1:size(source,3)
                    tz(j,:)=ir.register(double(source(:,:,j)));
                end
            end
            cumulative_t{center_z_index+direction*dz}...
                = cumulative_t{center_z_index+direction*(dz-1)} + tz(:,1:2);
            i = center_z_index+direction*dz;
            for ch = 1:numChannels
                aligned(:,:,ch,i) = BilinearImageRegistrator.shift(correctedZStack(:,:,ch,i),cumulative_t{i});
            end
        end
    end
    %%
    % StackViewer(aligned)
end