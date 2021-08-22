function [fn, ffn,sort_keys] = fastdir(dir_name,expression,opt,sort_key)
    if(nargin<2)
        expression = '';
    end
    if(nargin<3)
        opt = '';
    end
    if(nargin<4)
        sort_key = '';
    end
    jif = java.io.File(dir_name);
    if(~jif.isDirectory())
        warning('Not a directory');
        fn = {};
        ffn = {};
        return
    end
    jifl = jif.listFiles();
    fn = cell(jifl.length(),1);
    selected = false(jifl.length(),1);
    for i=1:jifl.length()
        fn{i} = char(jifl(i).getName());
        selected(i) = (isempty(expression) || ~isempty(regexp(fn{i},expression,'once','start')))...
                    &&(isempty(opt) || (strcmp(opt,'d') && jifl(i).isDirectory() ));
    end
    fn = fn(selected);
    switch(sort_key)
        case 'date'
            sort_keys = cellfun(@(x)java.io.File(fullfile(dir_name,x)).lastModified,fn);
        otherwise
            sort_keys = fn;
    end
    [~,i] = sort(sort_keys);
    fn = fn(i);
    sort_keys = sort_keys(i);
    
    if(nargout>1)
        ffn = cellfun(@(x)fullfile(dir_name,x),fn,'UniformOutput',false);
    end
end