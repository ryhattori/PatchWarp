function value = lininterp2_fast(V, x, y)
% linear interpolation, given set of a 2-d matrix and an x, y query
% This function acts as a reduced version of intep2 for bi-linear
% interpolation. Please note the algorithm used was retrieved from 
% https://en.wikipedia.org/wiki/Bilinear_interpolation -> unit square
% algorithm. 
%
% Input variables:
% V:            [mxn] [double] a 2-d array ( grayscale image ) 
% x:            [1xn] x ( column ) of location to be retrieved,
% y:            [1xn] y ( row ) of location to be retrieved
%
% Output:
% OUT:          [1xn] The intensity value at the x,y loc.
%--------------------------------------
% $ Ver: 1.3, 13/5/2012,  released by Georgios D. Evangelidis.
% Email: georgios.evangelidis@inria.fr
% $ Ver: 2.0, 12/5/2017,  released by Georgios D. Evangelidis.
% dan.erez.sd@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
    [y,x] = ind2sub(size(V),x);
end

%% prep 
x0 = floor(x);
x1 = ceil(x);
y0 = floor(y);
y1 = ceil(y);
szV = size(V);
x = mod(x,1);
y = mod(y,1);

%% find valid outputs 
validInds = x0> 0 & x1 < szV(2) & y0 > 0 & y1 < szV(1);
x(~validInds) = [];
y(~validInds) = [];
x0(~validInds) = [];
y0(~validInds) = [];
x1(~validInds) = [];
y1(~validInds) = [];

%% calc near by values
f00 = V(sub2ind(szV,y0,x0));
f01 = V(sub2ind(szV,y1,x0));
f10 = V(sub2ind(szV,y0,x1));
f11 = V(sub2ind(szV,y1,x1));

%% calc 
valueOfValidPts = f00.*(1-mod(x,1)).*(1-y)+f10.*x.*(1-y)+f01.*(1-x).*y+f11.*x.*y;

%% deal with poitns out of range
value = nan(numel(validInds),1);
value(validInds) = valueOfValidPts;




end