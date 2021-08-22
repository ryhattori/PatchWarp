function out = spatial_interp_patchwarp(in, warp, transform, nx, ny)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUT = SPATIAL_INTERP(IN, WARP, STR, TRANSFORM, NX, NY)
% This function implements the 2D spatial interpolation of image IN 
%(inverse warping). The coordinates defined by NX,NY are projected through 
% WARP thus resulting in new subpixel coordinates. The intensity values in 
% new pixel coordinates are computed via bilinear interpolation
% of image IN. 
%
% Input variables:
% IN:           the input image which must be warped,
% WARP:         the warp transform,
% TRANSFORM:    the type of adopted transform: {'translation','euclidean','affine','homography'}
% NX:           the x-coordinate values of horizontal side of ROI (i.e. [xmin:xmax]),
% NY:           the y-coordinate values of vertical side of ROI (i.e. [ymin:ymax]),
%
% Output:
% OUT:          The warped (interpolated) image
%--------------------------------------
% $ Ver: 1.3, 13/5/2012,  released by Georgios D. Evangelidis.
% Email: georgios.evangelidis@inria.fr
% $ Ver: 2.0, 12/5/2017,  released by Georgios D. Evangelidis.
% dan.erez.sd@gmail.com
% 
% Modified by Ryoma Hattori for PtachWarp
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

[xx, yy] = meshgrid(nx, ny);
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
out = lininterp2_fast(in, xy_prime(1,:)', xy_prime(2,:)');
% if fast_warping == 1
%     out = lininterp2_fast(in, xy_prime(1,:)', xy_prime(2,:)');
% else
%     out = interp2(in, xy_prime(1,:), xy_prime(2,:), str);
% end

out(isnan(out))=0;%replace Nan
out=reshape(out,length(ny),length(nx));
