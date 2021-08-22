function [Y, d] = max2d_subpixel(Z)
% Assume the curvature around peak is symmetric!
%
%% Check input
assert(isequal(size(Z),[3,3]),'Wrong matrix size. Input matrix should be numerical 3x3 type.');
if(any(isnan([Z(1,2) Z(2,1), Z(2,2),Z(2,3),Z(3,2)])))
    warning('Cannot handle NaN');
    Y = Z(2,2);
    d = [0 0];
    return;
end

assert(max(Z([2 4 5 6 8]))==Z(5),'Wrong matrix. Should be centered around its maximum.');

% pinv([2 -1 -1 1;1 -1 0 1; 2 -1 1 1; 1 0 -1 1; 0 0 0 1; 1 0 1 1; 2 1 -1 1; 1 1 0 1; 2 1 1 1])
% pinv([1 -1 0 1; 1 0 -1 1; 0 0 0 1; 1 0 1 1; 1 1 0 1;])

p = [1/4 1/4 -1 1/4 1/4; -1/2 0 0 0 1/2;0 -1/2 0 1/2 0;0 0 1 0 0 ] ...
    *[Z(1,2);Z(2,1);Z(2,2);Z(2,3);Z(3,2)];

if(p(1)>=0)
    Y=Z(2,2);
    d = [0 0];
    return
end

d = [-p(2)/(2*p(1)) -p(3)/(2*p(1))];
Y = p(4) - ( p(2)^2+p(3)^2) / (4*p(1));

% [y,x] = meshgrid(-2:0.1:2,-2:0.1:2);
% z = p(1)*(x.^2+y.^2)+p(2)*x+p(3)*y+p(4);
% imagesc(z);

end