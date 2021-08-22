function imn = imnormalize(im,r1,r2,offset,c)
    f = fspecial('disk',r1);
    f2 = fspecial('disk',r2);
    if(isa(im,'gpuArray'))
        imsm = zeros(size(im),'gpuArray');
        imsm2 = zeros(size(im),'gpuArray');
    else
        imsm = zeros(size(im));
        imsm2 = zeros(size(im));
    end
    r2=ceil(r2);
    n_cells = conv2(ones(size(im,1),size(im,2)),f2,'same');
    for i=1:size(im,3)
        imsm(:,:,i) = conv2(double(im(:,:,i)),f,'same');
        imsm2(:,:,i) = conv2(double(im(:,:,i)),f2,'same')./n_cells;
    end
    imn = mean(im(:))*(imsm+offset)./(imsm2+offset);
    if(nargin>4 && c)
        imn_new = imn;
        imn_new(:) = 0;
        [x, y]=meshgrid(1:size(imn,1),1:size(imn,2));
        for i=1:size(imn,3)
            [centers,radii] = imfindcircles(imn(:,:,i),(r1+r2)/2,'Sensitivity',0.98);
            for j=1:size(centers,1)
                imn_new(((x-centers(j,1)).^2 + (y-centers(j,2)).^2) < radii(j)^2)=1;
            end
        end
        imn = imn_new;
    end
end