%VISUALIZE Aplies envelope detection and log-compression (if desired) and
%plots the image. All settings are defined in the visParamHandler object
%obj.visParams (see >> doc visParamHandler).
%
%INPUT:
% - imgData: compounded high-frequency image dataset on the grid defined 
%            by z and x (no envelope detection, no log-compression). Can be
%            RF or IQ.
% - Z: 1D array or 2D grid with z-positions of the image pixels [m]
% - X: 1D array or 2D grid with x-positions of the image pixels [m]
%
%OUTPUT:
% - imgData: compounded image on the grid defined by Z and X after envelope 
%            detection and log-compression (if not deactivated). 
%
%AUTHOR: 
%hans-martin.schwab@web.de (May2019)

function imgData = visualize(obj, imgData, Z, X)
    %display
    %   Z,X optional
    
% start computation time:
ticVis = tic;
    
% rename params for shorter code:
vp = obj.visParams;
    
% get data size:
Nz = size(imgData,1);

%% Data Processing:

% compound images if uncompounded:
imgData = sum(imgData,3);

if vp.applyEnvelope
    % apply envelope:
    if isreal(imgData)
        % get absolute of analytic signal for real data:
        imgData = abs(RF2IQ_BL(imgData, round(Nz*obj.visParams.padFactHilbert)));
        % crop to original size that dataset had before padding:
        imgData = imgData(1:Nz,:);
    else
        % get absolute for IQ data:
        imgData = abs(imgData);
    end
end

% norm to maximum of entire data set (do not norm if maximum is 0!):
maxVal = max(abs(imgData(:)));
if maxVal~=0
    imgData = imgData./maxVal;
end

% apply log-compression only if .dynamicRange is not set to []:
if ~isempty(vp.dynamicRange)

    % if no envelope detection was applied, apply neg and pos logcompression:
    if vp.applyEnvelope
        % apply log compression:
        imgData = 20*log10(imgData);
        % apply dynamic range (crop off below and above thresholds):
        imgData(imgData<-vp.dynamicRange(1)) =  -vp.dynamicRange(1);
        if ~isscalar(vp.dynamicRange)
            imgData(imgData>-vp.dynamicRange(2)) =  -vp.dynamicRange(2);
        end
    else
        % force real (might be IQ):
        imgData = real(imgData);
        % mask negative values:
        negValMask = (imgData<0);
        % posify all values:
        imgData(negValMask) = -imgData(negValMask);
        % apply log compression and shift to positiv values:
        imgData = 20*log10(imgData) + vp.dynamicRange(1);
        % apply dynamic range value cropping:
        imgData(imgData<0) =  0;
        % negify negativ values:
        imgData(negValMask) = -imgData(negValMask);
    end
    
end

%% Dipsplay: 

% plot figure if figure handle is provided:
if checkFigHandle(vp.figHandle)
    
    % use different plot function depending on grid type:
    if vp.isEquidistGrid
        if nargin<3
            % display data on equidistant carthesian grid of pixels:
            imagesc(real(imgData)); axis normal; 
            xlabel('x [pixels]'); ylabel('z [pixels]');  
            
        elseif isvector(Z)
            % display data on equidistant carthesian grid:
            imagesc(X*1e3, Z*1e3, real(imgData)); axis image; 
            xlabel('x [mm]'); ylabel('z [mm]');   
            
        else
            % display data on equidistant carthesian grid:
            imagesc(X(1,:)*1e3, Z(:,1)*1e3, real(imgData)); axis image;                
            xlabel('x [mm]'); ylabel('z [mm]');   
            
        end
    else
        % display data on non-equidistant grid:
        pcolor(X*1e3, Z*1e3, real(imgData)); set(gca,'Color', 'k'); shading flat; 
        axis image ij; xlabel('x [mm]'); ylabel('z [mm]'); 
            
    end
    
    % apply further plot properties:
    colormap gray; colorbar; 
    if ~isempty(vp.figTitle), title(vp.figTitle); end
    
    % draw figure immediately (if not supressed)
    if ~vp.supressDrawnow, drawnow; end 
    
end
   
% stop computation time:
obj.compTimeVis = toc(ticVis);

% print computation time if desired:
if obj.printCompTime
    fprintf('computation time "visualization": %.3fs \n',obj.compTimeVis);
end

end  

%% LOCAL FUNCTIONS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert RF to IQ (analytic) signal with bandlimitation to avoid artifacts:
function data = RF2IQ_BL(data, N)

    % get size before FT:
    Ndata = size(data,1);

    % transfrom into frequency domain: 
    data = fft(data, N, 1);

    % crop lower sideband to avoid redundant computations:
    data = data(1:ceil(N/2), :, :);

    % weight data at zero-frequency to get analytic signal after IFT of upper
    % sideband:
    data(1,:,:) = 0;%.5*data(1,:,:);

    % apply band limitation filter as multiplicatoin in FD:
    data = bsxfun(@times, data, tukeywin(ceil(N/2),.1));

    % apply IFT and scale by 2 to get analytic signal:
    data = 2*ifft(data, Ndata, 1);

end
