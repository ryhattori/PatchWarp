function aligned = process_zstack_an_RH(fileNames)

addpath(genpath('/oasis/tscc/scratch/ryhattori/aki_rank_MC_RH'))
% addpath('Z:\People\Assaf\shield\MC__ZforRyoma')
A=process_zstack_RH(fileNames);
mip = max(A, [], 4);

majorpath=pwd;

% mkdir max_projection
% cd max_projection

B=squeeze(A);
stackfilename=[fileNames(1:end-7) 'stack.tif'];
write_tiff(stackfilename,int16(B(:,:,:)));


sumfilename=[fileNames(1:end-7) 'sum.tif'];
imwrite(uint16(mip),sumfilename)

cd(majorpath)


end