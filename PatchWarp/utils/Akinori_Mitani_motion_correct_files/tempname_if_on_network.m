function temp_fn = tempname_if_on_network(fn)
    persistent isNetworkDrive
    temp_fn = '';
    if(isunix && ~ismac)
        temp_fn = tempname('/dev/shm/');
    elseif(ispc)
        if(isempty(isNetworkDrive))
            isNetworkDrive = containers.Map('KeyType','char','ValueType','logical');
            drives = java.io.File('').listRoots();
            for i=1:numel(drives)
                isNetworkDrive(char(drives(i))) = ...
                   strcmp('Network Drive',char(javax.swing.filechooser.FileSystemView.getFileSystemView().getSystemTypeDescription(drives(i)))); 
            end
        end
        ffn = char(java.io.File(fn).getAbsoluteFile());
        if(length(ffn)>2 && isNetworkDrive.isKey(ffn(1:3)) && isNetworkDrive(ffn(1:3)))
            temp_fn = tempname;
        end
    end
end