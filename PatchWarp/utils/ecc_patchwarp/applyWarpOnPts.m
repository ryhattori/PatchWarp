function [valueAtInds] = applyWarpOnPts(inds,imIn,warp,transform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%valueAtInds = SPATIAL_INTERP(inds,imIn,warp,transform)
% This function find the intensity value of points after they have been 
% deplaced by a warp/transformation. 
%
% Input variables:
% inds:         the input image which must be warped,
% imIn:         the image to be sampled,
% WARP:         the warp transform,
% TRANSFORM:    the type of adopted transform: {'translation','euclidean','affine','homography'}
%
% Output:
% OUT:          The intensity value at the warped locations of inds.
%--------------------------------------
% $ Ver: 1.3, 13/5/2012,  released by Georgios D. Evangelidis.
% Email: georgios.evangelidis@inria.fr
% $ Ver: 2.0, 12/5/2017,  released by Georgios D. Evangelidis.
% dan.erez.sd@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% in affine or euclidean case, make the warp matrix 3x3
if (strcmp(transform,'affine')||strcmp(transform,'euclidean'))
   if size(warp,1)==2
       warp=[warp;zeros(1,3)];
   end
end

if strcmp(transform,'translation')
    warp = [eye(2) warp];
    warp = [warp; zeros(1,3)];
end

[yy,xx] = ind2sub(size(imIn),inds);
xy=[xx(:)';yy(:)';ones(1,length(yy(:)))];

%3x3 matrix transformation
A = warp;
A(3,3) = 1;

% new coordinates
xy_prime = A * xy;

if strcmp(transform,'homography')

    % division due to homogeneous coordinates
    xy_prime(1,:) = xy_prime(1,:)./xy_prime(3,:);
    xy_prime(2,:) = xy_prime(2,:)./xy_prime(3,:);
end

% Ignore third row
xy_prime = xy_prime(1:2,:);

% Subpixel interpolation
valueAtInds = lininterp2_fast(imIn, xy_prime(1,:)', xy_prime(2,:)');

valueAtInds(isnan(valueAtInds))=0;%replace Nan

end