function [template, selected] = make_template(stack,reduction_ratio,max_iter,threshold)
    if(nargin<2 || isempty(reduction_ratio))
        reduction_ratio = [16 16];
    end
    if(numel(reduction_ratio)<2)
        reduction_ratio = [reduction_ratio reduction_ratio];
    end
    
    if(nargin<3||isempty(max_iter))
        max_iter = 5;
    end
    
    if(nargin<4||isempty(threshold))
        threshold = 0.2;
    end
    
    m = size(stack,1);
    n = size(stack,2);
    N = size(stack,3);
    
    sm = ceil(m/reduction_ratio(1));
    sn = ceil(n/reduction_ratio(2));

    small_stack = zeros(sm,sn,N);
    for i = 1:N
        small_stack(:,:,i) = imresize(stack(:,:,i),[sm sn]);
    end
    
    data_2d = reshape(small_stack,[sm*sn N]);
    selected = true(N,1);
    for iter = 1:max_iter
        mean_2d = mean(data_2d(:,selected),2);
        
        c = zeros(N,1);
        for i = 1:N
            r = corrcoef(data_2d(:,i),mean_2d);
            c(i) = r(1,2);
        end
        
        selected = c>quantile(c,1-threshold);
    end
    
    template = mean(stack(:,:,selected),3);
end
