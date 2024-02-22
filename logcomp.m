% function imgDB = logcomp(imgHF,DR,padFactor)
%
% applies an envelope detection and a logarithmic compression on the image
% "imgHF", whereas the dynmic range of the output image "imgDB" is determined
% by "DR" in dB (optional, default:40). If DR is a two-element vector, the
% first element will be treated as lower and the second as upper threshold.
% The dimensions of "imgHF" and "imgDB" are in the order (axial,lateral,frame).
% "padFactor" determines the factor with which the data set is padded for
% the fft in the envelope detection to avoi wrapping (optional, default:0).
%
% author: Hans-Martin Schwab

function imgDB = logcomp(imgHF,DR,padFactor)

% set default for dynmaic range:
if nargin<2
   DR = [40,0]; 
end
if isscalar(DR), DR = [DR,0]; end

% set default for pad factor:
if nargin<3
   padFactor = 1; 
end

% apply envelope (for complex data, IQ data is assumed):
if isreal(imgHF)
    imgDB = abs(RF2IQ_BL(imgHF ,size(imgHF,1)*round(padFactor) ));
else
    imgDB = abs(imgHF);
end

% crop to original size before padding:
imgDB = imgDB(1:size(imgHF,1),:,:);
                
% norm to maximum of entire data set:
imgDB = imgDB./max(imgDB(:));

% apply log compression:
imgDB = 20*log10(imgDB);

% apply dynamic range:
imgDB(imgDB<-DR(1)) =  -DR(1);
imgDB(imgDB>-DR(2)) =  -DR(2);

end

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
    data = data .* tukeywin(ceil(N/2),.1);

    % apply IFT and scale by 2 to get analytic signal:
    data = 2*ifft(data, N, 1);
    data = data(1:Ndata,:,:);
end