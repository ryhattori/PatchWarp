function d = downsample_mean(x,n,dim)
    if(nargin<3)
        dim = [];
    end
    d = downsample_chunk(x,n,dim,'mean');
end
