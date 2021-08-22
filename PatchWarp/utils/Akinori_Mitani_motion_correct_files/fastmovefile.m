function s = fastmovefile(source,dist)
    dist_parent = fileparts(dist);
    if(~java.io.File(dist_parent).exists)
        mkdir(dist_parent);
    end
    s = java.io.File(source).renameTo(java.io.File(dist));
end