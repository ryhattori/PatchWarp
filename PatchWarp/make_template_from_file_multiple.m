function [template,selected] = make_template_from_file_multiple(fn, stack_range, align_ch, n_ch, threshold, fftdenoise, to_save, network_temp_copy)
    if(nargin<1)
        fn = [];
    end
    if(nargin<2)
        align_ch = [];
    end
    if(nargin<3)
        threshold = [];
    end
    if(nargin<4)
        to_save = true;
    end
   
    
    if(isempty(fn))
        [filename, pathname] = uigetfile('*.tif','Select a file to open');
        if(isequal(0,filename) || isequal(0,pathname))
            template = [];
            return
        else
            fn = fullfile(pathname,filename);
        end
    end
    
    
    [pathstr,filename,~] = fileparts(fn{stack_range(1),:});
    dirname = fullfile(pathstr,'template');
    if(~java.io.File(dirname).isDirectory())
        mkdir(dirname);
    end
%     filename = filename(1:end-4);
%     for i = 1:length(stack_range)
%         filename = [filename, '_', num2str(stack_range(i))];
%     end
    
    fn_template_tif = fullfile(dirname,[filename '_avg.tif']);
    fn_template_mat = fullfile(dirname,[filename '_avg.mat']);
    
    info = cell(length(stack_range),1);
    stacks = cell(length(stack_range),1);
    for i = 1:length(stack_range)
         [stacks{i},info{i}, ~] = read_tiff(fn{stack_range(i)}, align_ch, n_ch, network_temp_copy);
%        [stacks(i),~,~,info{i}] = read_SI4_image(fn{stack_range(i)},chs);
    end
    
    if(any(cellfun(@isempty,stacks)))
        error('Specified channel was not in the file.');
    end
    if(isempty(stacks))
       error('isempty'); 
    end
    
    stacks_class = class(stacks{1});
    stacks = cat(3, stacks{:});
%     [stacks,~,~,info] = read_SI4_image(fn,chs);
   
    
    template = zeros(size(stacks,1),size(stacks,2), 1,stacks_class);
    selected = cell(1,1);
    for i = 1:1
        [template(:,:,i),selected{i}] = make_template_fftdenoise(stacks,[],[],threshold, fftdenoise);
    end
    
    if(to_save)
        write_tiff(fn_template_tif,template,info);
        save(fn_template_mat,'template','selected','fn');
    end
end