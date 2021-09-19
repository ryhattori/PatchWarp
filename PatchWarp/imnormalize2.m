function imn = imnormalize2(im,r)
    f = fspecial('disk',r);
    if(isa(im,'gpuArray'))
        imsm = zeros(size(im),'gpuArray');
    else
        imsm = zeros(size(im));
    end
    n_cells = conv2(ones(size(im,1),size(im,2)),f,'same');
    for i=1:size(im,3)
        imsm(:,:,i) = conv2(double(im(:,:,i)),f,'same')./n_cells;
    end
    imn = mean(im(:))*im./imsm;
end