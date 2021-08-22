function [template, t] = make_stable_average_RH(stack,align_ch,use_mex)
    if(~iscell(stack))
        stack = {stack};
    end
    if(nargin<2)
        align_ch = numel(stack);
    end

    target = int16(mean(stack{align_ch},3));
    source = int16(stack{align_ch});
    
    if use_mex
        t_tmp = mexBilinearRegistrator(target,source);
    else
        ir = BilinearPyramidImageRegistrator(target,0.75,3);
        t_tmp = zeros(size(source,3),2);
        for j = 1:size(source,3)
            t_tmp(j,:)=ir.register(double(source(:,:,j)));
        end
    end

    
    
    t_tmp = t_tmp(:,1:2);
    t0 = zeros(1,2);
    for j=1:2
        tmp = sort(t_tmp(:,j));
        window = ceil(0.2*size(t_tmp,1));
        [~,k] = min(tmp((window+1):end)-tmp(1:(end-window)));
        t0(j) = (tmp(window+k)+tmp(k))/2;
    end
    t = bsxfun(@minus,t_tmp,t0);
    template = zeros(size(stack{1},1),size(stack{1},2),numel(stack));
    for ii = 1:numel(stack)
        source = int16(stack{ii});
        if use_mex
            shifted = mexBilinearShift(source,t);
        else
            shifted = zeros(size(source),'single');
            for j = 1:size(source,3)
                shifted(:,:,j)=BilinearPyramidImageRegistrator.shift(source(:,:,j),t(j,:));
            end
        end
        template(:,:,ii) = mean(shifted,3);
    end
end

%%
