%FILTER Applies the filter processing steps to the receive data, which can
%contain a bandfilter, conversion from RF to IQ and upsampling if desired.
%All settings are defined in the filtParamHandler object obj.filtParams 
%(see >>doc filtParamHandler).
%
%INPUT:
% - RXData: 2D or 3D array with receive data with temporal samples
%           of all elements in columns and different transmission
%           events in the third dimension
%
%OUTPUT:
% - RXData: 2D or 3D array with receive data with temporal samples
%           of all elements in columns and different transmission
%           events in the third dimension after filtering, IQ conversion
%           and upsampling
% - dtu:  temporal spacing of RXData [s]. The temporal spacing might not 
%         equal 1/fs, because the data might have been upsampled to 
%         increase the acuracy of DAS. This information has to be
%         transferred to .reconstruct by this value.
%
%AUTHOR: 
%hans-martin.schwab@web.de (May2019)
       
function [RXData, dtu] = filter(obj, RXData)
%filter
%   filterer optional

% start computation time:    
ticFilt = tic;

%% get parameters from input:

% get filter params for shorter code:
fp = obj.filtParams;

% get size:
Nt = size(RXData,1);

% get sample number after upsampling:
Ntu = ceil(Nt*fp.upsampleFactor);

% return temporal spacing of upsampled data:
dtu = 1/obj.fs/obj.filtParams.upsampleFactor;

%% compute filter coeffs:
   
% initiliaze filtCoeffsFD as empty:
filtCoeffsFD = [];

% compute filter coefficients according to filter type:
switch fp.filtType
    
    % apply no filter (and hence no padding):
    case 'off'
        Nt_pad = Nt;
        filtCoeffsFD = ones(Nt_pad,1);
        
     % apply hanning filter with -6dB bandwidth defined in .param:
    case 'simpleBP'
        assert(~isempty(obj.fc),['The property .fc of the class DASReconstructor '...
               'must be defined to use the filter option ''simpleBP''']);
        % get parameters:
        bw = fp.params;
        % apply fix lengths TD zeropadding (wrapping will still occur for strong band limitation):
        Nt_pad = Nt + 20;
        % get sample number of center frequeny:
        Ntc = obj.fc/obj.fs*Nt_pad;
        % get sample number of bandwidth:
        Nbw = Ntc*bw;
        % set filterCoeffs as hanning window:
        arg_angle = 2*pi*((0:ceil(Nt_pad/2)-1)-Ntc)/(2*Nbw); % get argument for cos
        arg_angle = min(max(arg_angle,-pi),pi); % apply hanning limits to argument
        filtCoeffsFD = .5*(1+cos(arg_angle)); % get hanning filter

    % apply a filter using fir1 that damps frequencies beyond the nyquist 
    % limit of the image spacing: 
    case 'antialias' 
        assert(~isempty(obj.fc),['The property .fc of the class DASReconstructor '...
               'must be defined to use the filter option ''simpleBP''']);
        % get parameters:        
        dr = fp.params(1);
        bwFact = fp.params(2); 
        Nt_filt = fp.params(3);  
        % compute aliasing filter in time domain:
        fc_frac = obj.fc/(obj.fs/2);
        bw_frac = 1/(2*dr/obj.c0) / obj.fs;
        bandLims_frac =  min(max(fc_frac + bwFact*bw_frac*[-.5,.5],.015),.95);
        filtCoeffsTD = fir1(Nt_filt, bandLims_frac);
    
    % apply a filter with manually set time-domain filter coefficients
    case 'manualTD' 
        filtCoeffsTD = fp.params;  
        
    % apply a filter with a manually set frequency-domain window function
    case 'manualFD' 
        % do not apply TD zeropadding (wrapping will occur):
        Nt_pad = Nt;
        % get FD filter coefficients:
        if isa(fp.params, 'function_handle')
            % get coefficients for param as function handle:
            filtCoeffsFD = fp.params(obj.fs/Nt_pad*(0:ceil(Nt_pad/2)-1));
        else
            % get coefficients for param as array by resampling:            
            filtCoeffsFD = interp1(0:length(fp.params)-1, fp.params,...
                                   linspace(0,length(fp.params)-1,ceil(Nt_pad/2)));
        end
end

% if FD coeffs were not defined yet, TD coeffs need to be transformed:
if isempty(filtCoeffsFD)
    % apply TD zeropadding (avoid wrapping):
    Nt_pad = Nt + length(filtCoeffsTD);
    % transform filter into FD:
    filtCoeffsFD = fft(filtCoeffsTD(:),Nt_pad,1);
    % crop lower sideband:
    filtCoeffsFD = filtCoeffsFD(1:ceil(Nt_pad/2));
    % shift phase to make filter zerophase:
    f = (0:ceil(Nt_pad/2)-1)'/Nt_pad;
    filtCoeffsFD = filtCoeffsFD.*exp(2*pi*1i*f*(length(filtCoeffsTD)-1)/2);% <- is that always the same as: filtCoeffsFD=abs(filtCoeffsFD); 
end

% get upsampled padded element number:
%(important for correct IFT with upsampling but without wrapping)
Ntu_pad = ceil(Nt_pad*fp.upsampleFactor);

%% plot debug figure:

% plot figure if figure handle is provided:
if checkFigHandle(fp.figHandle)
    
    % get normalized spectrum of RFData: 
    dspect = mean(mean(abs(fft(double(RXData),Nt_pad,1)),2),3); dspect = dspect/max(dspect);
    dspect = dspect(1:ceil(Nt_pad/2));
    
    % get approximized frequency vector (not accurate):
    f = linspace(0, obj.fs/1e6/2, ceil(Nt_pad/2));
    
    if isempty(filtCoeffsFD)
        % plot mean spectra:
        plot(f,dspect); xlabel('Frequency [MHz]'); title('Spectrum'); 
        legend('RXData');
    else
        % get normalized spectrum of filter: 
        fspect = abs(filtCoeffsFD); fspect = fspect/max(fspect);
        plot(f,dspect, f,fspect); xlabel('Frequency [MHz]'); 
        title('Data and filter spectra'); 
        legend('RXData','Filter');        
    end
    
    % plot band limits if provided:
    if 0%IF FC AND BW ARE GIVEN:
        hold all; 
        plot(max(fc_frac-.5*bw_frac,0)*[1,1]*fs/2e6,[0,1],'k--'); 
        plot(min(fc_frac+.5*bw_frac,1)*[1,1]*fs/2e6,[0,1],'k--'); 
        hold off;
    end
    
    % further plot settings:
    ylim([0,1]);

end

%% apply attenuation compensation:
% (compensation neglects transmit delays and frequency dependence)

% apply attenuation compensation only if not constant:
if obj.filtParams.attenCoeff ~= 0
    
    % get travel time equivalent two-way distance axis [cm]: 
    dist_tau = (0:Nt-1)' * obj.c0/obj.fs/2 * 1e2;    
    % get depth-amplification function:
    depthAmp = exp(dist_tau * log(10^(obj.filtParams.attenCoeff/20)));    
    % apply depth dependent amplification:
    RXData = bsxfun(@times, RXData, depthAmp); 
    
end

%% data processing in Fourier domain (upsampling, filtering, hilbert trafo):

% transfrom into frequency domain: 
RXData = fft(RXData, Nt_pad, 1);

% crop lower sideband to avoid redundant computations:
RXData = RXData(1:ceil(Nt_pad/2), :, :);

% weight data at zero-frequency to get analytic signal after IFT of upper
% sideband:
RXData(1,:,:) = 0;%.5*data(1,:,:);

% apply filter only if filter is not turned off:
if ~strcmp(fp.filtType,'off')
    % apply filter as multiplicatoin in FD:
    RXData = bsxfun(@times, RXData, filtCoeffsFD(:));
end

% apply IFT onto padded length and scale by 2 to get analytic signal:
%(also scaling by resample factor in case of upsampling required)
RXData = max(1,fp.upsampleFactor)*2*ifft(RXData, max(Ntu_pad,Nt_pad), 1);
% crop off padded part:
RXData = RXData(1:max(Ntu,Nt),:,:);
%NOTE: (FT interpolation is applied in case of upsampling (max(Ntu,Nt)=Ntu), 
% in case of downsampling (max(Ntu,Nt)=Nt), resampling will applied after 
% down mixing to avoid aliasing)
%TODO!!!

% get real part of analytic signal if 'RF' data is desired:
if strcmp(fp.outputFormat,'RF')
    RXData = real(RXData);
end

%% apply upmixing if signal is demodulated:

% % apply downmixing if desired:
% if strcmp(fp.inputFormat,'IQDemod');
%     
%     % create spacing and time axis:
%     dt = 1/sampleFreq_RF/max(upsampleFactor,1);
%     t = dt*(0:(max(Ntu,Nt)-1))';
%     % (spacing and number of samples depends on whether it is upsampling
%     %  ,which has already happened, or downsampling which will happen
%     %  after downmixing to avoid aliasing)
%     
%     % apply parallelized downmix operation:
%     RXData = bsxfun(@times, RXData, exp(-1i*2*pi*t*centerFreq));
% 
% end

%% resampling:

% apply downsampling if upsampleFactor<0:
NskipSmpl = floor(max(1/fp.upsampleFactor,1));
RXData = RXData(1:NskipSmpl:end,:,:);

% stop computation time:
obj.compTimeFilt = toc(ticFilt);

% print computation time if desired:
if obj.printCompTime
    fprintf('computation time "filter":        %.3fs \n',obj.compTimeFilt);
end

end