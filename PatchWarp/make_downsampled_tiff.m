function processed = make_downsampled_tiff(save_dir,opt)
    if(nargin<2)
        opt = '';
    end
    
    fns = dir(fullfile(save_dir,'*_summary.mat'));
%	fns = dir(fullfile(save_dir,'*.tif'));
    if(isempty(fns))
        fprintf('no summary files in the folder %s\n',save_dir);
        downsampled = [];
%        return;
    end
    downsampled_c = cell(1,1,length(fns));
    downsampled_perstack_c = cell(1,1,length(fns));
    for i=1:length(fns)
        tmp=load(fullfile(save_dir,fns(i).name));
        downsampled_c{i} = tmp.downsampled;
        if(isfield(tmp,'downsampled_perstack'))
            downsampled_perstack_c{i} = tmp.downsampled_perstack;
        end
    end
    downsampled = cell2mat(downsampled_c);
    downsampled_perstack = cell2mat(downsampled_perstack_c);

    ffn = fullfile(save_dir,'downsampled',sprintf('downsampled_%d.tif',tmp.n_downsampled));
    ffn_max = fullfile(save_dir,'downsampled',sprintf('downsampled_%d_max.tif',tmp.n_downsampled));
    ffn_mean = fullfile(save_dir,'downsampled',sprintf('downsampled_%d_mean.tif',tmp.n_downsampled));
    ffn_align = fullfile(save_dir,'downsampled','downsampled_perstack.tif');
    
    processed = exist(ffn,'file') & exist(ffn_align,'file');
    if(strcmp(opt,'-s'))
        return;
    end
    
    warning off;
    mkdir(fullfile(save_dir,'downsampled'));
    warning on;
    
    try
        write_tiff(ffn,int16(downsampled));
        write_tiff(ffn_max,int16(max(downsampled, [], 3)));
        write_tiff(ffn_mean,int16(mean(downsampled, 3)));
    catch
    end
    
    if(~isempty(downsampled_perstack))
       write_tiff(ffn_align,int16(downsampled_perstack));
    end

end
