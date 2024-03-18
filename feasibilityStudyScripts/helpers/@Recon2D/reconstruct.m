%RECONSTRUCT Applies all the reconsruction processing steps and returns a
%compounded HF image. All settings are defined in the reconParamHandler object
%obj.reconParams (see >> doc reconParamHandler).
%
%INPUT:
% - RXData: 2D or 3D array with receive data with temporal samples
%           of all elements in columns and different transmission
%           events in the third dimension, which will be compounded.
% - Z: 1D array or 2D grid with z-positions of the image pixels [m]
% - X: 1D array or 2D grid with x-positions of the image pixels [m]
% - dt: temporal spacing of RXData [s]. If .filter was applied prior to this
%       method, the temporal spacing might not equal 1/fs, because the 
%       data might have been upsampled to increase the acuracy of DAS. 
%OUTPUT:
% - imgData: compounded high-frequency image dataset on the grid defined 
%            by z and x (no envelope detection, no log-compression). Can be
%            RF or IQ, depending on if "RXData"  is RF or IQ.
%AUTHOR: 
%hans-martin.schwab@web.de (May2019)
       
function imgData = reconstruct(obj, RXData, Z, X, dt)
    %reconstruct
    %   Z,X optional?
    
% start computation time:    
ticRecon = tic;

% short name for recon param handler for denser code:
rp = obj.reconParams;
    
%% input checking:

% check if sizes of .TXPulseDelays and .elemPos are consistent with RXData:
obj.checkConsistency(RXData);

% get input data size:
[Nt,Na_RX,Nacq] = size(RXData);

% parse axes inputs:
if isvector(Z) && isvector(X)
   % keep axes as vectors only if arrayfun will be executed on GPU:
   if isa(RXData,'gpuArray') || isa(Z,'gpuArray')
       % transpose vectors into the respective dimension if currently wrong:
       if size(Z,1)==1, Z = Z(:); end
       if size(X,2)==1, X = X(:)'; end
   else
       % convert axes into grids if arrayfun will be executed on the CPU 
       % (would return an error otherwise bc no singular dimension expansion):
       [Z,X] = ndgrid(Z,X);
   end
else
    % check consistency of grid size:
    assert(all(size(X)==size(Z)), ['The input arguments ''x'' ''z'' must',...
           'be vectors or full grids of equal size!']) 
end
% turn axes inot GPU arrays if RXData is a gpu array
if isa(RXData,'gpuArray') && ~isa(Z,'gpuArray')
    Z = gpuArray(Z);
    X = gpuArray(X);
end

%% Prepare general settings:

% convert element position to non-property array (faster):
elemPos = obj.elemPos;

% if defined use TX aperture, otherwise TXAperture equals RX aperture:
elemPosTX = obj.elemPosTX;
if isempty(elemPosTX), elemPosTX = elemPos; end
Na_TX = size(elemPosTX,1);

% convert pulse delays to distance:
elemPulseDelays = obj.TXPulseDelays*obj.c0;

% determine if acquisition is PA (photoacsoutics) or US (ultrasound)
% (for PA, only one row of TX delays is provided, for US one per element):
isPAMode = size(elemPulseDelays,1)==1;

% converter to samples:
dist2smpl = 1/(dt*obj.c0);

%% Prepare interpolation settings:

% translate interpolation method from string to index (1,2 or 3):
interpMethod = find(strcmp(rp.interpMethod, {'next','linear','phaseNext','phaseLinear'}));

% precompute constant phase factor in exponential shift term for phase interp:
if ~isreal(RXData) % only if IQ
    fullPhaseShift = -1i*2*pi * obj.fc*dt; 
else
    % if RXData is not IQ data, there must be no phase shift:
    fullPhaseShift = 0;
end

%% Prepare apodization settings:

% get apodization settings: 
maxAngle = rp.apodCutoffAngle;
apodWindow = rp.apodWeights;
N_win = length(apodWindow);

% set booleans to decide on apodization options in recon:
applyApodWeights = ~isempty(rp.apodWeights);
applyHanningWeights = strcmp(rp.apodWeights,'hanning');
if applyHanningWeights,apodWindow=0; end % no chars as nested vars for GPU
%NOTE: even if the variables "apodWindow" and "N_win" are not reached in 
% recon, they need to be initialized, otherwiese the GPU returns an error! 

% check if number of element smeets size of normal vectors:
if length(obj.normalAngleRX) == Na_RX  
    % apply normal angles from object:
    normalAngleRX = obj.normalAngleRX;
elseif isscalar(obj.normalAngleRX)
    % get all element normal angles if only scalar is given:
    normalAngleRX = obj.normalAngleRX*ones(Na_RX,1);    
else
    error('The size of ".normalAngleRX" does not match number of columns in "RXData"');
end

% convert angles to tangents (to avoid computations of tan in kernel): 
lowerLimTan = tand(-rp.apodCutoffAngle);
upperLimTan = tand( rp.apodCutoffAngle); 

% get sin and cosines of negative normal angles to rotate pixel-element
% angle in realitve orientation to normal angle in kernel to avoid 
% problems with unambigious tangent function:  
sinTh = sind(-normalAngleRX);
cosTh = cosd(-normalAngleRX);

%% Prepare insonification weight settings:
% Note that in PA (photoacoustics) mode, there is not TX insonification. Hence,
% the weighting is disabled and the parameter ".insoniWeightType" is ignored!

% set insonification weight settings:
if strcmp(rp.insoniWeightType,'off') || isPAMode
    
        % do not apply mask calculation in reconstruction by setting 
        % TXDelayCutoff to a false boolean (will neglect computation of mask):
        TXDelayCutoff      = false;
        applyInsoniWeights = false; 
        
    elseif strcmp(rp.insoniWeightType,'WLBinary')    
        
        % check if .fc exists:
        assert(~isempty(obj.fc),['The property .fc of the class DASReconstructor '...
         'must be defined to use the insonification weight option ''WLBinary''']);  
     
        % get cutoff value as centerwavelength times cutoff fraction in [m]:
        TXDelayCutoff      = rp.insoniWeightParams(1)*obj.c0/obj.fc;
        applyInsoniWeights = false;
                
   elseif strcmp(rp.insoniWeightType,'WLConti') % NAME MACHT KEINEN SINN!  
       
        % check if .fc exists:       
        assert(~isempty(obj.fc),['The property .fc of the class DASReconstructor '...
         'must be defined to use the insonification weight option ''WLConti''']); 
     
        % get cutoff value as centerwavelength times cutoff fraction in [m]:
        TXDelayCutoff      = rp.insoniWeightParams(1)*obj.c0/obj.fc;
        applyInsoniWeights = true;
        
    else
        error('insonification weight type %s  not supported!',rp.insoniWeightType);
end    

%% apply recon:

% apply TX event weighting (synthetic TX apodization) if desired:
if ~isempty(obj.reconParams.acqWeights)
    
   % create hanning window if no array but a string 'hanning' is provided:
   if strcmp(obj.reconParams.acqWeights,'hanning')
        acqWeights = hanning(obj.Nacq);
   else
        acqWeights = rp.acqWeights;
        assert(numel(obj.reconParams.acqWeights)==obj.Nacq,...
               ['The property ".reconParams.acqWeights" must be empty or have a',...
                'number of elements that matches the third dimaension of "RXData"']);        
   end
   % apply one weight to each frame in RX dataset:
   RXData = bsxfun(@times, RXData, permute(acqWeights(:),[3,2,1])); 
   
end

% for GPU-timeout issues, there is an option to process each frame
% individually (in series):
if rp.applySerialFrameComp || ~rp.applyCompounding
    
   % for serial frame computation, imgData is initilaized and filled up:
   if rp.applyCompounding
       % for compounding only one frame is initialized:
       imgData = zeros(size(Z,1),size(X,2)) * RXData(1);% (*X(1) to get correct data type)
   else
       % for non-compounding one frame per acquisition is initialized:
       imgData = zeros(size(Z,1),size(X,2),obj.Nacq)*RXData(1);    
   end
   
   for i = 1:Nacq
    % indizes which are iterated over are set to the same frame:
    %(iAcq__start and iAcq__end must be defined individually, probably due 
    % to a bug in MATLAB gpu handling. One array with indizes does not work)     
    iAcq_start = i;
    iAcq_end   = i;
    % compound acquisitions is desired:
    if rp.applyCompounding
        % compute reconstruction as nested function in paralellized function call:     
        imgData = imgData + arrayfun(@singlePixelDAS, Z, X); 
    else
        % compute reconstruction as nested function in paralellized function call:  
        imgData(:,:,i) = arrayfun(@singlePixelDAS, Z, X); 
    end
   end
   
else
    
   % include all frames in one call: 
   %(iAcq__start and iAcq__end must be defined individually, probably due 
   % to a bug in MATLAB gpu handling. One array with indizes does not work)
   iAcq_start = 1;
   iAcq_end   = Nacq;
   % compute reconstruction as nested function in paralellized function call:
   imgData = arrayfun(@singlePixelDAS, Z, X);
   
end

%% plot grid:

if checkFigHandle(rp.figHandle)
    
    % plot grid:
    hold off;
    if isvector(X), [Z,X] = ndgrid(Z,X); end
    scatter(X(:)*1e3,Z(:)*1e3,'.','MarkerEdgeColor',[.5,.5,.5]);
    legendText = {'pixel locations'};
    hold all;
    
    % plot transducer:
    scatter(obj.elemPos(:,2)*1e3,obj.elemPos(:,1)*1e3,'g.');  
    legendText{end+1} = 'transducer elements';
    
    % plot TX transducer (if different):
    if ~isempty(obj.elemPosTX)
        scatter(obj.elemPosTX(:,2)*1e3,obj.elemPosTX(:,1)*1e3,'m.');
        legendText{end}   = 'transducer elements (RX)';
        legendText{end+1} = 'transducer elements (TX)';
    end

    if isscalar(obj.normalAngleRX)
       theta = obj.normalAngleRX*ones(obj.Na,1);
    else
       theta = obj.normalAngleRX; 
    end
    
    % set length of the line plot as 5x element distance:
    lineLength = 4*sqrt(sum(diff(obj.elemPos(1:2,:),1,1).^2))*1e3;
    % plot normal angles of elements:
    plot(obj.elemPos(:,2)'*1e3 + lineLength*[zeros(obj.Na,1),sind(theta)]',...
         obj.elemPos(:,1)'*1e3 + lineLength*[zeros(obj.Na,1),cosd(theta)]','g');       
    
    hold off;
    % apply further plot properties:
    xlabel('x [mm]'); ylabel('z [mm]'); colormap gray; axis image ij;
    title('Reconstruction Grid'); legend(legendText,'Location','south');
    drawnow; 
    
end

%% timing:

% stop computation time:
obj.compTimeRecon = toc(ticRecon);

% print computation time if desired:
if obj.printCompTime
    fprintf('computation time "reconstruct":   %.3fs \n',obj.compTimeRecon);
end

%% NESTED FUNCTIONS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general single pixel reconstruction (GPU compatible): 
function pixVal = singlePixelDAS(zPos, xPos)

    % initialize value:
    pixVal = 0*RXData(1);
    w_apod = 0;
    w_interp = 0; % must be initialized for GPU, even if unused!
    
    % initialize sumThisPixel as true (might be set false if mask is used):
    addThisPixel = true;
    
    % Loop over all frames:
    for i_acq = iAcq_start:iAcq_end
%     for i_acq = 1:Nacq
     
        % initialize TXDelay and TXDelay_i:
        TXDelay = inf;
        TXDelay_i = inf;
        TXDelayDiff = inf;

        % compute TX delays depending on US mode or PA (photoacoustics) mode: 
        if isPAMode
        
            % in PA mode, the TX delay equals the laser pulse delay
            % stored in "elemPulseDelays" as single row entry:
            TXDelay = elemPulseDelays(1,i_acq);
                    
        else
            
            % in US mode, the TX delays are determined depending on element
            % distance and transmit delays of the elements:
            
            % sum over active TX aperture:
            for i_a = 1:Na_TX

                % compute TX delay:
                elemTXDelay = sqrt((elemPosTX(i_a,1)-zPos).^2 ...
                                        +(elemPosTX(i_a,2)-xPos).^2);

                % determine both TX delay and delay difference if TXDelayCutoff
                % was not set to false (if false no insonification mask is applied):
                if TXDelayCutoff  

                    % get delay of last element as current delay before
                    % recomputation of current delay:
                    TXDelay_iMinus1 = TXDelay_i;

                    % update curretn TX Delay including pulse delay:
                    TXDelay_i = elemPulseDelays(i_a,i_acq) + elemTXDelay;                

                    % set new delay and new difference to current if delay is
                    % smaller than the currently smallest delay:
                    if TXDelay_i < TXDelay
                        TXDelay = TXDelay_i;
                        TXDelayDiff = TXDelay_i - TXDelay_iMinus1;
                    end       

                else

                    % get minimum TX delay including pulse Delay:
                    TXDelay = min(TXDelay, elemPulseDelays(i_a,i_acq) + elemTXDelay);

                end

            end

            % correction and mask computation if TXDelayCutoff is not false:
            if TXDelayCutoff

                % if TXDelay is delay of 1st element, the computation of TOF
                % difference to the previous element was false, the inf-init
                % was preserved and the difference to next element is used:
                % ! MUST BE ADAPTED TO COVER inf in TXPulseDelays ! 
                if isinf(TXDelayDiff)
                    TXDelayDiff = sqrt((elemPosTX(1,1)-zPos).^2+(elemPosTX(1,2)-xPos).^2)...
                                 +elemPulseDelays(1,i_acq)-elemPulseDelays(2,i_acq)...
                                 -sqrt((elemPosTX(2,1)-zPos).^2+(elemPosTX(2,2)-xPos).^2);
                end

                % determine if pixel is insonified enough to be included in mask: 
                addThisPixel = abs(TXDelayDiff) < TXDelayCutoff;

            end
            
        end
        
        % apply binary insonification mask (neglect out of mask pixels):
        if addThisPixel
            
            % sum over RX aperture:
            for i_a = 1:Na_RX 
                    
                % get hypothenuses of distances between pixel and element:
                diffZ = zPos-elemPos(i_a,1);
                diffX = xPos-elemPos(i_a,2);
                
                % get scalar product of elem-pixel vector and normal vector
                % to check if pixel is behind aperture:
                cosRotated = -diffX*sinTh(i_a)+diffZ*cosTh(i_a);
                
                % get relative tangent of angle by rotating hypothenuses 
                % by element normal angle to compare to cutoff tangent:
                elemTan = (diffX*cosTh(i_a)+diffZ*sinTh(i_a))/cosRotated;                    

                % apply angle dependended binary apodization: 
                if cosRotated>0 && elemTan>lowerLimTan && elemTan<upperLimTan
                    
                    % compute RX delay:
                    elemRXDelay = sqrt(diffZ.^2 + diffX.^2);
                    
                    % compute two-way delay in samples and then as integer NN index:
                    DASDelaySmpl = (TXDelay + elemRXDelay)*dist2smpl; % +1 because indexes in Matlab are 1-based not 0-based!
                    DASDelayIdx  = ceil(DASDelaySmpl);
                    
                    % only add to sum if value within RFData:
                    if DASDelayIdx<Nt && DASDelayIdx>1
                        
                        % compute interpolation weight for non next neighbor
                        % interpolation:
                        if interpMethod > 1, w_interp = DASDelayIdx - DASDelaySmpl; end
                        
                        % add pixel val without weights or with if existing:
                        if applyApodWeights 
                            if applyHanningWeights
                                % apodisation Weight: (2*pi/(2*maxAngle)=pi/maxAngle)
                                w_apod = .5+.5*cosd(180/maxAngle*atand(elemTan));
                            else
                                % get apodization weigth index:
                                i_w = ceil(N_win * (elemTan-lowerLimTan)/...
                                           (upperLimTan-lowerLimTan));
                                % compute apodization weight:
                                w_apod = apodWindow(i_w);
                            end
                            % add apodization weighted value to sum:
                            if interpMethod==1 % NEXT NEIGHBOR
                                % sum apodization weighted values:
                                pixVal =  pixVal + w_apod*RXData(DASDelayIdx, i_a, i_acq); 
                            elseif interpMethod==2 % LINEAR
                                % sum apodization weighted and interpoolation weighted values:
                                pixVal = pixVal + w_apod...
                                        *( (w_interp-1)*RXData(DASDelayIdx  , i_a, i_acq)...
                                          + w_interp   *RXData(DASDelayIdx-1, i_a, i_acq));
                            elseif interpMethod==3 % PHASESHIFT (FROM NEXT NEIGHBOR)
                                % sum apodization weighted and interpolation shifted values: 
                                pixVal =  pixVal + w_apod*RXData(DASDelayIdx, i_a, i_acq)...
                                                  *exp(fullPhaseShift*w_interp);
                            elseif interpMethod == 4 % PHASESHIFT (LINEAR FROM 2 NEIGhBORS)
                                % sum interpolation shifted values: 
                                pixVal =  pixVal + w_apod...
                                         *( (1-w_interp)*RXData(DASDelayIdx  , i_a, i_acq)*exp(fullPhaseShift*w_interp)...
                                           + w_interp   *RXData(DASDelayIdx-1, i_a, i_acq)*exp(fullPhaseShift*(w_interp-1)));                            
                            end        
                        else
                            % add apodization weighted value to sum (w_apod=1 would be slow):
                            if interpMethod == 1 % NEXT NEIGHBOR
                                % sum values:
                                pixVal =  pixVal + RXData(DASDelayIdx, i_a, i_acq); 
                            elseif interpMethod == 2 % LINEAR
                                % sum  interpoolation weighted values:
                                pixVal = pixVal + ...
                                        +(1-w_interp)*RXData(DASDelayIdx  , i_a, i_acq)...
                                        + w_interp   *RXData(DASDelayIdx-1, i_a, i_acq);
                            elseif interpMethod == 3 % PHASESHIFT (FROM NEXT NEIGHBOR)
                                % sum interpolation shifted values: 
                                pixVal =  pixVal + RXData(DASDelayIdx, i_a, i_acq)...
                                                  *exp(fullPhaseShift*w_interp);
                            elseif interpMethod == 4 % PHASESHIFT (LINEAR FROM 2 NEIGhBORS)
                                % sum interpolation shifted values: 
                                pixVal =  pixVal + ...
                                         +(1-w_interp)*RXData(DASDelayIdx  , i_a, i_acq)*exp(fullPhaseShift*w_interp)...
                                         + w_interp   *RXData(DASDelayIdx-1, i_a, i_acq)*exp(fullPhaseShift*(w_interp-1));                            
                            end                      
                        end

                    end

                end

            end
        
            % apply insonification weight if pixel within mask and mask not binary:  
            if applyInsoniWeights
                % compute continuous pixel weight for insonification weighting:
                pixVal = pixVal*(TXDelayCutoff-abs(TXDelayDiff))/TXDelayCutoff;
            end
            
        end

    end
end
%%Note: It would be faster to compute the apodization cutoffs once and then
%%iterate over Nacq for the sum, but this is uncompatible with arrayfun, 
%%because then the TX delays had to be computed for each RX element because
%%no TX delay arrays can be stored. 
    
end