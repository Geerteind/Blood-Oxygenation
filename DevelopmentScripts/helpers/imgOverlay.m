% function [ paus_img ] = imgOverlay( us_data, pa_data, threshold, compression, colmap )
% Creates an Pa/Us overlay image
%
% OUTPUT:
%
%   paus_img: rgb image with PA image overlayed on US image
%
% INPUT:
%
%   us_data:    2-D matrix, dataset containing reconstructed ultraound 
%               image data (can be log compressed or linear).
%
%   pa_data:    2-D matrix, dataset containing reconstructed photoaocustic 
%               image data (can be log compressed or linear).
%
%   threshold:  scalar between 0 and 1, fractianal intensity value of PA 
%               image below which PA image is completely tansperent.
%
%   compression:scalar between 0 and inf, determines the relation between 
%               PA image intensity and transparency for values higher than 
%               threshold. If set to 0, no transparency is used. If set to 
%               1, transparency increases linearly with image intensity. 
%
%   colmap:     Nx3 matrix, colormap for PA overlay image (default: hot). 
%               The US image is grayscale by default. If you want to 
%               change the US colormap, you can set colmap as a cell array
%               colmap = {colmap_pa,colmap_us}.
%
% Author: Hans-Martin Schwab (hans-martin.schwab@web.de / hans-martin.schwab@rub.de)



function [ paus_img ] = imgOverlay( us_data, pa_data, threshold, compression, colmap )

% check input arguments:
if nargin<3
    threshold=.5;
    compression=1;
    pa_cmap=hot(64);
    us_cmap=gray(64);
elseif nargin<4
    compression=1;
    pa_cmap=hot(64);
    us_cmap=gray(64);
elseif nargin<5
    pa_cmap=hot(64);
    us_cmap=gray(64);
else
    if iscell(colmap)
        pa_cmap=colmap{1}; 
        us_cmap=colmap{2};
    else
        pa_cmap=colmap; 
        us_cmap=gray(size(pa_cmap,1));
    end
end
n_indicees_pa_cmap=size(pa_cmap,1);
n_indicees_us_cmap=size(pa_cmap,1);

% If PA and US do not have the same size, US size is used:
if any(size(us_data)~=size(pa_data))
   pa_data=imresize(pa_data,size(us_data)); 
end

% pa_cmap=flip(autumn(64),1);
% pa_cmap=autumn(64);
% pa_cmap=hot(64);
% pa_cmap=cool(64);
% pa_cmap=hot(64); pa_cmap(50:end,:)=[];

% norm data and apply threshold:
us_data=us_data-min(us_data(:));
us_data=us_data/max(us_data(:));
pa_data=pa_data-min(pa_data(:));
pa_data=pa_data/max(pa_data(:));
    pa_data(pa_data<threshold)=threshold;
    pa_data=pa_data-threshold;
    pa_data=pa_data/max(pa_data(:));
%     pa_data=pa_data.^compression;

% convert grayscale to indexed image
[pa_ind, pa_map] = gray2ind(pa_data,n_indicees_pa_cmap);
[us_ind, map] = gray2ind(us_data,n_indicees_us_cmap);

% convert indexed to rgb image:
pa_img=ind2rgb(pa_ind,pa_cmap);
us_img=ind2rgb(us_ind,us_cmap);%

% compute transparency:
alpha=pa_data;
% alpha(alpha<=threshold)=threshold;
% alpha=alpha-threshold;
% alpha=alpha/max(alpha(:));
alpha=alpha.^compression;

% compute rgb overlay:
alpha=repmat(alpha,[1,1,3]);
paus_img=us_img.*(1-alpha)+pa_img.*alpha;

% imshow(pa_img)
% imshow(us_img)
% imshow(paus_img)


end