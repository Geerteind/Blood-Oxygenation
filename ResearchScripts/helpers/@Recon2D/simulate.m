%SIMULATE performs a time-of-flight based simulation for the given
%geometry and a set of point-scatterers and returns an RF dataset.
%For PA acquisitions (.TXPulseDelays has only one row), the
%point scatterers are treated as point sources.
%This method is faster than .simulateSerial, but the timing is less
%accurate, because delays are only extact to 1/.fs/16.
%INPUTS:
% - scatPos: (Nscat x 2) or (Nscat x 3) matrix containing the
%            positions of Nscat scatterers [z,x] with amplitude 1 or
%            positions and amplitudes [z,x,amp] in rows
% - nPwr: normalized scalar power of additive Gaussian noise on the 
%         output RF data (relative to maximum TX signal power 1)
% - Nt: scalar positive number of samples in output RF dataset
% - upsmplFct: (optional) scalar positive factor by which the temporal 
%              signal is upsampled before TOFs are rounded. Increasing
%              this value increases the precision of TOF timing (default:16)
%OUTPUT:
% - RXData: (Nt x Na x Nacq) array with simulated RF data assuming a 
%           Gauss-pulse transmission from the TX element positions
%           according to the transmit delays, single scattering at 
%           the scatterpositions and receive from the RX elements.
%
%TODO: 
% - inputchecking and defaults for inputs
% - GET RID OF LOOP OVER CHANNELS! 


function RXData = simulate(obj, scatPos, nPwr, Nt, upsmplFct)
   
    narginchk(4,5);
    %input check: <...>
    
    % parse inputs:
    scatPos = gather(scatPos); % sim does not work with gpu arrays
    if isempty(obj.elemPosTX),TXPos=obj.elemPos; else,TXPos=obj.elemPosTX; end 
    if nargin<5, upsmplFct=16; end
    
    % initialize RXData:
    RXData = zeros(Nt, obj.Na, obj.Nacq);
    
    % initialize upsampled excitation pulse:
    t0 = gauspuls('cutoff', obj.fc);
    excitPulse = gauspuls(-t0 : 1/obj.fs/upsmplFct : t0, obj.fc);
    
    % init scatter intensity (does not change if undefined):
    scatIntensity = 1;
    scatIntensity_i = 1;

    % get scatterer intensity:
    if size(scatPos,2)==3, scatIntensity = scatPos(:,3); end
        
    % init permuted scatPos:
    scPos = permute(scatPos(:,1:2),[3,2,1]);
    
    % check if object is in PA (photoacoustics) mode:
    isPAmode = size(obj.TXPulseDelays,1)==1;
    
    % get delays from scatterer to receive elements:
    RXScatDelays = sqrt(sum((obj.elemPos-scPos).^2,2))/obj.c0;               
    RXScatDelays = permute(RXScatDelays,[3,1,2]);
        
    %iterate over all transmission events:
    for i_tx = 1:obj.Nacq
        
        % get delay when point emits wave depending on PA or US mode:
        if isPAmode
            % get laser pulse delay defined in .TXPulseDelays:
            TXScatDelay = obj.TXPulseDelays(1, i_tx);
        else
            % get delay when wavefront reaches scatterer:
            TXScatDelay = min(obj.TXPulseDelays(:,i_tx) + ...
              sqrt(sum((TXPos-scPos).^2,2))/obj.c0,[],1);
        end
        
        %*** insert pulse amplitude on each channel at respective delay:   
        for i_chan = 1:obj.Na
            
            % get closest temporal index to center of received pulse:
            i_t = round((TXScatDelay(:)+RXScatDelays(:,i_chan)) * obj.fs * upsmplFct);
            
            % delete indices and intensities that are out of range:
            if size(scatPos,2)==3, scatIntensity_i = scatIntensity(i_t>0 & i_t<=Nt*upsmplFct); end
            i_t = i_t(i_t>0 & i_t<=Nt*upsmplFct);
            
            % create dataset with upsampled data to fill with intensities:
            rawLineData = zeros(Nt*upsmplFct,1);
            rawLineData(i_t) = rawLineData(i_t) + scatIntensity_i;
            
            % convolve rawLineData with exictation pulse:
            rawLineData = conv(rawLineData, excitPulse, 'same');
            
            % add noise and fill into RXData array:
            RXData(:,i_chan,i_tx) = rawLineData(1:upsmplFct:end);
            
        end
            
    end

    % add noise to RX data if defined:
    if nPwr~=0, RXData = RXData + sqrt(nPwr)*randn(Nt, obj.Na, obj.Nacq); end

end


%% Alternative code for parallel channel computation:
% (turned out to be slower)

%         %*** insert pulse amplitude on all channels at respective delay:
%
%         % get closest temporal index to center of received pulse:
%         i_t = round((TXScatDelay(:)+RXScatDelays) * obj.fs * upsmplFct);
% 
%         % delete indices and intensities that are out of range:
%         if size(scatPos,2)==3
%             scatIntensity_i = scatIntensity * ones(1,obj.Na);
%             scatIntensity_i = scatIntensity_i(~(i_t<=1 & i_t>Nt*upsmplFct)); 
%         end
%         i_t(i_t<=1 & i_t>Nt*upsmplFct) = nan;
%         i_t = i_t + Nt*upsmplFct*(0:obj.Na-1);
%         i_t = i_t(~isnan(i_t));
%         
%         % create dataset with upsampled data to fill with intensities:
%         rawRXData = zeros(Nt*upsmplFct*obj.Na,1);
%         rawRXData(i_t) = rawRXData(i_t) + scatIntensity_i;
% 
%         % convolve rawLineData with exictation pulse:
%         rawRXData = imfilter(reshape(rawRXData,[Nt*upsmplFct,obj.Na]), excitPulse(:));
% 
%         % add noise and fill into RXData array:
%         RXData(:,:,i_tx) = rawRXData(1:upsmplFct:end,:);
