function [results, warp, warpedImage, rho_final] = ecc_patchwarp(image, template, levels, noi, transform, delta_p_init, learningRate, ptsPerIter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ECC image alignment algorithm
%[RESULTS, WARP, WARPEDiMAGE] = ECC(IMAGE, TEMPLATE, LEVELS, NOI, TRANSFORM, DELTA_P_INIT)
%
% This m-file implements the ECC image alignment algorithm as it is
% presented in the paper "G.D.Evangelidis, E.Z.Psarakis, Parametric Image Alignment
% using Enhanced Correlation Coefficient.IEEE Trans. on PAMI, vol.30, no.10, 2008"
%
% ------------------
% Input variables:
% IMAGE:        the profile needs to be warped in order to be similar to TEMPLATE,
% TEMPLATE:     the profile needs to be reached,
% NOI:          the number of iterations per level; the algorithm is executed
%               (NOI-1) times
% LEVELS:       the number of levels in pyramid scheme (set LEVELS=1 for a
%               non pyramid implementation), the level-index 1
%               corresponds to the highest (original) image resolution
% TRANSFORM:    the image transformation {'translation', 'euclidean', 'affine', 'homography'}
% DELTA_P_INIT: the initial transformation matrix for original images (optional); The identity
%               transformation is the default value (see 'transform initialization'
%               subroutine in the code). In case of affine or euclidean transform,
%               DELTA_P_INIT must be a 2x3 matrix, in homography case it must be a 3x3 matrix,
%               while with translation transform it must be a 2x1 vector.
% learningRate: The rate at which each iteration's effect on the warp is
%               applied ( 0 - 1 = 0% - 100% ). default = .75. If your registration 
%               often diverges reduce this downwards (  .1 - .9 are typical values ). 
% ptsPerIter:   The number of points sampled in each iteration - default
%               150. It is best to have at least 15 - 25 points per paramenter
%               ( i.e. homogrpahy - > 8 parms -> 150 points ).  
%
% For example, to initialize the warp with a rotation by x radians, DELTA_P_INIT must
% be [cos(x) sin(x) 0 ; -sin(x) cos(x) 0] for affinity
% [cos(x) sin(x) 0 ; -sin(x) cos(x) 0 ; 0 0 1] for homography.
%
%
% Output:
%
% RESULTS:   A struct of size LEVELS x NOI with the following fields:
%
% RESULTS(m,n).warp:     the warp needs to be applied to IMAGE at n-th iteration of m-th level,
% RESULTS(m,n).rho:      the enhanced correlation coefficient value at n-th iteration of m-th level,
% WARP :              the final estimated transformation [usually also stored in RESULTS(1,NOI).WARP ].
% WARPEDiMAGE:        the final warped image (it should be similar to TEMPLATE).
%
% The first stored .warp and .rho values are due to the initialization. In
% case of poor final alignment results check the warp initialization
% and/or the overlap of the images.
% -------------------
% $ Ver: 1.4, 12/2/2013,  released by Georgios D. Evangelidis, INRIA, FRANCE
% For any comment, please contact the author
% Email: georgios.evangelidis@inria.fr, evagelid@ceid.upatras.gr
% $ Ver: 2.0, 12/5/2017,  released by Georgios D. Evangelidis.
% dan.erez.sd@gmail.com
%
% This software is provided "as is" without any kind of warranty. Also, it
% is provided for research purposes only. In any case, please cite the above paper.
% 
% Modified by Ryoma Hattori for PtachWarp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


break_flag=0;

if nargin<5
    error('-> Not enough input arguments');
end

if ~exist('learningRate','var')
    learningRate = .75;
end
if ~exist('ptsPerIter','var')
    ptsPerIter = 8*15;
end



transform = lower(transform);
if ~(strcmp(transform,'affine')||strcmp(transform,'euclidean')||strcmp(transform,'homography')||strcmp(transform,'translation'))
    error('-> Not a valid transform string')
end

sZi3 = size(image,3);
sZt3 = size(template,3);
initImage = image;
initTemplate = template;
if sZi3>1
    if ((sZi3==2) || (sZi3>3))
        error('Unknown color image format: check the number of channels');
    else
        image=rgb2gray(uint8(image));
    end
end

if sZt3>1
    if ((sZt3==2) || (sZt3>3))
        error('Unknown color image format: check the number of channels');
    else
        template = rgb2gray(uint8(template));
    end
end

template = double(template);
image = double(image);

%% pyramid images
% The following for-loop creates pyramid images in cells IM and TEMP with varying names
% The variables IM{1} and TEMP{1} are the images with the highest resoltuion

% Smoothing of original images
% f = fspecial('gaussian',[7 7],.5);
f = get7x7Gaussfilt();
TEMP{1} = conv2(template,f,'same');
IM{1} = conv2(image,f,'same');
TEMP{1} = template;
IM{1} = image;

% create pyramid 
for nol=2:levels
    IMCurr = conv2(IM{nol-1},ones(2)./4);
    IM{nol} = IMCurr(1:2:end,1:2:end);
    TEMPCurr = conv2(TEMP{nol-1},ones(2)./4);
    TEMP{nol} = TEMPCurr(1:2:end,1:2:end);
end


%% transform initialization

% In case of translation transform the initialiation matrix is of size 2x1:
%  delta_p_init = [p1;
%                  p2]
% In case of affine transform the initialiation matrix is of size 2x3:
%
%  delta_p_init = [p1, p3, p5;
%                  p2, p4, p6]
%
% In case of euclidean transform the initialiation matrix is of size 2x3:
%
%  delta_p_init = [p1, p3, p5;
%                  p2, p4, p6]
%
% where p1=cos(theta), p2 = sin(theta), p3 = -p2, p4 =p1
%
% In case of homography transform the initialiation matrix is of size 3x3:
%  delta_p_init = [p1, p4, p7;
%                 p2, p5, p8;
%                 p3, p6,  1]
if strcmp(transform,'translation')
    nop=2; %number of parameteres
    if nargin==5;
        warp=zeros(2,1);
    else
        if (size(delta_p_init,1)~=2)|(size(delta_p_init,2)~=1)
            error('-> In translation case the size of initialization matrix must be 2x1 ([deltaX;deltaY])');
        else
            warp=delta_p_init;
        end
    end
end
if strcmp(transform,'euclidean')
    nop=3; %number of parameteres
    if nargin==5;
        warp=[1 0 0; 0 1 0; 0 0 0];
    else
        if (size(delta_p_init,1)~=2)||(size(delta_p_init,2)~=3)
            error('-> In euclidean case the size of initialization matrix must be 2x3');
        else
            warp=[delta_p_init;zeros(1,3)];
        end
    end
end
if strcmp(transform,'affine')
    nop=6; %number of parameters
    if nargin==5;
        warp=[1 0 0; 0 1 0; 0 0 0];
    else
        if (size(delta_p_init,1)~=2)|(size(delta_p_init,2)~=3)
            error('-> In affine case the size of initialization matrix must be 2x3');
        else
            warp=[delta_p_init;zeros(1,3)];
        end
    end
end

if strcmp(transform,'homography')
    nop=8; %number of parameteres
    if nargin==5;
        warp=eye(3);
    else
        if (size(delta_p_init,1)~=3)|(size(delta_p_init,2)~=3)
            error('-> In homography case the size of initialization matrix must be 3x3');
        else
            warp=delta_p_init;
            if warp(3,3)~=1
                error('The ninth element of homography must be equal to 1');
            end
        end
    end
end

% in case of pyramid implementation, the initial transformation must be
% appropriately modified
for ii=1:levels-1
    warp=next_level(warp, transform, 0);
end

%% Run ECC algorithm for each level of pyramid
if numel(noi) < nol
    noi = repmat(noi,nol,1);
end
for nol=levels:-1:1
    
    im = IM{nol};
%     [vx1,vy1]=gradient(im);
    vx = conv2(double(im),[.5,0,-.5],'same');
    vy = conv2(double(im),[.5,0,-.5]','same');

    
    
    temp = TEMP{nol};
    
    [A,B]=size(temp);
    % Warning for tiny images
    if prod([A,B])<400
        disp(' -> ECC Warning: The size of images in high pyramid levels is quite small and it may cause errors.');
        disp(' -> To avoid such errors you could try fewer levels or larger images.');
        disp(' -> Press any key to continue.')
        pause
    end
    
    % Define the rectangular Region of Interest (ROI) by nx and ny (you can modify
    % the ROI). Here we just ignore some image margins. Margin is equal to
    % 5 percent of the mean of [height,width].
    m0=mean([A,B]);
    margin=floor(m0*.05/(2^(nol-1)));
    margin =0;%no-margin - modify these two lines if you want to exclude a margin
    nx=margin+1:B-margin;
    ny=margin+1:A-margin;
    temp=double(temp(ny,nx,:));
    
    %%%%   temp=temp-mean(temp(:)); % zero-mean image; is useful for brightness change compensation, otherwise you can comment this line
    % MODIFIED 12/1/2013 (we subtract the zero-mean inside the loop by taking
    % into account only the overlapping area
    
    %% check which parts of the image is even relevant 
    highGrad = find(vx > mean(vx(:)) | vy > mean(vy(:)));
    
    
    %% ECC, Forwards Additive Algorithm -------------------------------
    for i=1:noi(nol);

        %ADDED 12/05/2017; MODIFIED 12/05/2017
        % select random pts to apply the method onto
        indsToTest = randi([1,numel(highGrad)],ptsPerIter,1);
        indsToTest = highGrad(indsToTest) ;
        
        [indsToTestI,indsToTestJ] = ind2sub(size(im),indsToTest);
        temp_curr = temp(indsToTest);
        
        
        wim = applyWarpOnPts(indsToTest,im,warp,transform);
        
        
        
        %ADDED 7/8/2012; MODIFIED 12/2/2013
        % define a mask to deal with warping outside the image borders
        % (they may have negative values due to the subtraction of the mean value)
        %         ones_map = spatial_interp(ones(size(im)), warp, 'nearest', transform, nx, ny); %inverse (backward) warping
        ones_map = wim > 0 ;
        numOfElem = sum(sum(ones_map~=0));
        
        meanOfWim = sum(sum(wim.*(ones_map~=0)))/numOfElem;
        meanOfTemp = sum(sum(temp_curr.*(ones_map~=0)))/numOfElem;
        
        
        wim = wim-meanOfWim;% zero-mean image; is useful for brightness change compensation, otherwise you can comment this line
        tempzm = temp_curr-meanOfTemp; % zero-mean template
        
        wim(ones_map==0) = 0; % for pixels outside the overlapping area
        tempzm(ones_map==0)=0;
        
        
        %Save current transform
        if (strcmp(transform,'affine')||strcmp(transform,'euclidean'))
            results(nol,i).warp = warp(1:2,:);
        else
            results(nol,i).warp = warp;
        end
        
        results(nol,i).rho = temp_curr(:)'*wim(:) / norm(tempzm(:)) / norm(wim(:));
        
        if (i == noi) % the algorithm is executed (noi-1) times
            break;
        end
        
        % Gradient Image interpolation (warped gradients)
        wvx = applyWarpOnPts(indsToTest,vx,warp,transform);
        wvy = applyWarpOnPts(indsToTest,vy,warp,transform);
        
        
        % Compute the jacobian of warp transform
        J = warp_jacobian_onPts(indsToTestJ,indsToTestI, warp, transform);
        
        
        % Compute the jacobian of warped image wrt parameters (matrix G in the paper)
        G = image_jacobian_onPts(wvx, wvy, J, nop);
        
        
        
        % Compute Hessian and its inverse
        C= G' * G;% C: Hessian matrix
        con=cond(C);
        if con>1.0e+15
%             disp('->ECC Warning: Badly conditioned Hessian matrix. Check the initialization or the overlap of images.')
        end
        i_C = inv(C);
        
        % Compute projections of images into G
        Gt = G' * tempzm(:);
        Gw = G' * wim(:);
        
        
        %% ECC closed form solution
        
        % Compute lambda parameter
        num = (norm(wim(:))^2 - Gw' * i_C * Gw);
        den = (tempzm(:)'*wim(:) - Gt' * i_C * Gw);
        lambda = num / den;
        
        % Compute error vector
        imerror = lambda * tempzm - wim;
        
        % Compute the projection of error vector into Jacobian G
        Ge = G' * imerror(:);
        
        % Compute the optimum parameter correction vector

        delta_p = i_C * Ge * learningRate;
        
        
        if (sum(isnan(delta_p)))>0 %Hessian is close to singular
%             disp([' -> Algorithms stopped at ' num2str(i) '-th iteration of ' num2str(nol) '-th level due to bad condition of Hessian matrix.']);
%             disp([' -> Current results have been saved at results(' num2str(nol) ',' num2str(i) ').warp and results(' num2str(nol) ',' num2str(i) ').rho.']);
%             disp([' -> If you enabled a multilevel running, the output variables (warp, warpedImage) have been computed after mapping the current warp into the high-resolution level']);
            
            break_flag=1;
            break;
        end
        
        % Update parmaters
        warp = param_update(warp, delta_p, transform);
        
        
    end
    
    if break_flag==1
        break;
    end
    
    % modify the parameteres appropriately for next pyramid level
    if (nol>1)&(break_flag==0)
        warp = next_level(warp, transform,1);
    end
    
end


if break_flag==1 % this conditional part is only executed when algorithm stops due to Hessian singularity
    for jj=1:nol-1
        warp = next_level(warp, transform,1);
        %m0=2*m0;
    end
    %margin=floor(m0*.05);
    %nx=margin+1:size(template,2)-margin;
    %ny=margin+1:size(template,1)-margin;
end

if break_flag == 1
    final_warp = warp;
else
    final_warp = results(1,noi(1)).warp;
end

% return the final warped image using the whole support area (include
% margins)

nx2 = 1:B;
ny2 = 1:A;

for ii = 1:sZi3
    warpedImage(:,:,ii) = spatial_interp_patchwarp(double(initImage(:,:,ii)), final_warp, transform, nx2, ny2);
end

warpedImage = uint8(warpedImage);
warp = final_warp;

rho_final = results(end,end).rho;
%%%