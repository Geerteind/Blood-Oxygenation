%SIMULATEFULLTX please check .simulate for reference. This function does te
% same, but is more accurate and slower, because all TX waves are
% considered, not only the first one to arrive at the scatterer. 
% Use this simulation for focussed waves or other wavefronts that rely on
% destructive interference and if the sidelobes induced by TX play an
% essential role!
%
%AUTHOR: hans-martin.schwab@web.de (Jan21)

function RXData = simulateFullTX(obj, scatPos, nPwr, Nt, upsmplFct)
   
        narginchk(4,5);
        %input check: <...>
        scatPos = gather(scatPos); % sim does not work with gpu arrays
        if isempty(obj.elemPosTX),TXPos=obj.elemPos; else,TXPos=obj.elemPosTX; end 
        if nargin<5, upsmplFct=16; end
            
        % initialize RXData as zeros or noisy dataset:
        if nPwr == 0
            RXData = zeros(Nt, obj.Na, obj.Nacq);
        else
            RXData = sqrt(nPwr)*randn(Nt, obj.Na, obj.Nacq);            
        end
        
        % initialize upsampled excitation pulse:
        t0 = gauspuls('cutoff', obj.fc);
        excitPulse = gauspuls(-t0 : 1/obj.fs/upsmplFct : t0, obj.fc);
        
%         % get half window size (numberof samples within 3 x center wavelength):
%         dt = 1/obj.fs; 
%         Nt_win_half = floor(3 * 1/obj.fc / dt);
        % init scatter intensity (does ont change if undefined defined):
        scatIntensity = 1;
        scatIntensity_i = 1;
        
        % get scatterer intensity:
        if size(scatPos,2)==3, scatIntensity = scatPos(:,3); end
                        
        % init permuted scatPos:
        scPos = permute(scatPos(:,1:2),[3,2,1]);        
        
        % check if object is in PA (photoacoustics) mode:
        isPAmode = size(obj.TXPulseDelays)==1;

        % get delays from scatterer to receive elements:
        RXScatDelays = sqrt(sum((obj.elemPos-scPos).^2,2))/obj.c0;               
        RXScatDelays = permute(RXScatDelays,[3,1,2]);        
        
        %iterate over all transmission events:
        for i_tx = 1:obj.Nacq
            
            % get delay when point emits wave depending on PA or US mode:
            if isPAmode

                % get laser pulse delay defined in .TXPulseDelays:
                TXScatDelays = repmat(obj.TXPulseDelays(1, i_tx),obj.Na,1);

            else

                % get delay when wavefront reaches scatterer:
                TXScatDelays = obj.TXPulseDelays(:,i_tx) + ...
                  sqrt(sum((TXPos-scPos).^2,2))/obj.c0;
            end

            %*** insert pulse amplitude on each receive channel at respective delay:   
            for i_chan = 1:obj.Na
                
                % initialize empty raw line:
                rawLineData = zeros(Nt*upsmplFct,1);
                
                % loop over all TX elements:
                for i_a=1:obj.Na   

                    % get closest temporal index to center of received pulse:
                    i_t = round((permute(TXScatDelays(i_a,:,:),[3,1,2])...
                         +RXScatDelays(:,i_chan)) * obj.fs * upsmplFct);

                    % delete indices and intensities that are out of range:
                    if size(scatPos,2)==3, scatIntensity_i = scatIntensity(i_t>0 & i_t<=Nt*upsmplFct); end
                    i_t = i_t(i_t>0 & i_t<=Nt*upsmplFct);
                    
                    % fill dataset with upsampled data with intensities:
                    rawLineData(i_t) = rawLineData(i_t) + scatIntensity_i;
                    
                end

                % convolve rawLineData with exictation pulse:
                rawLineData = conv(rawLineData, excitPulse, 'same');

                % add noise and fill into RXData array:
                RXData(:,i_chan,i_tx) = RXData(:,i_chan,i_tx) + rawLineData(1:upsmplFct:end);

           end

        end
end


%                     % get time vectors for each receive element with 0 as receive time:
%                     t = bsxfun(@minus,(0:Nt-1)/obj.fs - TXScatDelays(i_a), RXScatDelays);
%                     % create noisy receive data as gauspuls around time zero for each element:
%                     RXData(:,:,i_tx) = RXData(:,:,i_tx) + scatIntensity*gauspuls(t, obj.fc)';
