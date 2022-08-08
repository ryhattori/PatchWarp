function [template, selected] = make_template_fftdenoise(stack,reduction_ratio,max_iter,threshold,fftdenoise)
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

    if (fftdenoise == true)
        freqmap_mean = fftshift(fft2(template));
        freqmap_mean_amp = log(abs(freqmap_mean));
        
        threshold = 0.5;
        freqmap_mean_amp_midlines0 = freqmap_mean_amp;
        freqmap_mean_amp_midlines0(round(m/2)-2:round(m/2)+2,:) = 0;
        power_max1 = max(freqmap_mean_amp_midlines0,[],1)';
        f_y = medfilt1(power_max1,10,'omitnan','truncate');
        power_max1_detrended = power_max1 - f_y;
%         [~,locs_max1] = findpeaks(power_max1_detrended,'MinPeakHeight', std_threshold*std(power_max1_detrended));
        locs_max1 = find(power_max1_detrended > threshold);
        
        freqmap_mean_amp_midlines0 = freqmap_mean_amp;
        freqmap_mean_amp_midlines0(:,round(n/2)-2:round(n/2)+2) = 0;
        power_max2 = max(freqmap_mean_amp_midlines0,[],2);
        f_y = medfilt1(power_max2,10,'omitnan','truncate');
        power_max2_detrended = power_max2 - f_y;
%         [~,locs_max2] = findpeaks(power_max2_detrended,'MinPeakHeight', std_threshold*std(power_max2_detrended));
        locs_max2 = find(power_max2_detrended > threshold);
        
        remove_mask = zeros(size(freqmap_mean_amp)); 
        remove_range = 3;
        if ~isempty(locs_max1)
            for i = 1:length(locs_max1)
                if (locs_max1(i)-remove_range) < 1
                    remove_mask(:, 1:locs_max1(i)+remove_range) = 1;
                elseif (locs_max1(i)+remove_range) > n
                    remove_mask(:, locs_max1(i)-remove_range:n) = 1;
                else
                    remove_mask(:, locs_max1(i)-remove_range:locs_max1(i)+remove_range) = 1;
                end
            end
        end
        if ~isempty(locs_max2)
            for i = 1:length(locs_max2)
                if (locs_max2(i)-remove_range) < 1
                    remove_mask(1:locs_max2(i)+remove_range, :) = 1;
                elseif (locs_max2(i)+remove_range) > m
                    remove_mask(locs_max2(i)-remove_range:m, :) = 1;
                else
                    remove_mask(locs_max2(i)-remove_range:locs_max2(i)+remove_range, :) = 1;
                end
            end
        end

        x=2;
        remove_mask(round(m/2)-x:round(m/2)+x, round(n/2)-x:round(n/2)+x) = 0;
        remove_mask = logical(remove_mask);
    
        temp = fftshift(fft2(template));
        temp(remove_mask) = 0;
        template = abs(ifft2(fftshift(temp)));
    end
end
