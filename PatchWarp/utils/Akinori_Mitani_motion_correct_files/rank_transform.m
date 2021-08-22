function result = rank_transform(im)
    result = zeros(size(im));
    for i=1:size(im,3)
        [~,~,ic] = unique(im(:,:,i));
        result(:,:,i) =  reshape(ic,size(im,1),size(im,2));
    end
end