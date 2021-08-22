function d = downsample_chunk(x,n,dim,method)
% d = downsample_chunk(x,n,dim,method)
   
    if(nargin<3 || isempty(dim))
        dim = 1;
    end
    if(nargin<4 || isempty(method))
        method = 'mean';
    end
    if(isinf(n))
        reshaped_x = x;
        ds = size(x);
        ds(dim)=1;
    else
        matdim = ndims(x);

        s = size(x);
        if(length(s)<dim)
            s(end+1:dim)=1;
            matdim=dim;
        end
        ds = s;
        ds(dim) = floor(s(dim)/n);

        trancated_x = subsref(x,struct('type',{'()'},'subs',{cat(2,repmat({':'},1,dim-1),{1:(n*ds(dim))},repmat({':'},1,matdim-dim))}));
        reshaped_x = reshape(trancated_x,[ds(1:(dim-1)) n ds(dim) ds((dim+1):end)]);
    end
    
    switch(method)
        case 'mean'
            d = reshape(mean(reshaped_x,dim),ds);
        case 'nanmean'
            d = reshape(nanmean(reshaped_x,dim),ds);
        case 'median'
            d = reshape(median(reshaped_x,dim),ds);
        case 'min'
            d = reshape(min(reshaped_x,[],dim),ds);
        case 'max'
            d = reshape(max(reshaped_x,[],dim),ds);
        otherwise
            error('unknown method');
    end
end
