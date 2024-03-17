%% GETUSIMAGE retrieves the ultrasound images from plane wave US receive data 
%
% INPUTS:
% - receiveDataRFUS: (Nt-by-Nel-by-Nfr) array that stores the plane 
%                    wave ultrasound  receive data for all frames
% - fs: scalar with temporal sample frequency of dataset [Hz]
% - x_elem: vector with trnasducer element positions in x-direction [m]
%           (z-positions are assumed 0)
% - c0: scalar with speed-of-sound used for reconstruction [m/s]
% - z: vector with z_positions of image pixels [m]
% - x: vector with x_positions of image pixels [m]
%
% OUPTUS:
% - imgUS: (Nz-by-Nx-by-Nfr) array with log compressed US images for each
%          frame scaled to maximum [dB]
%
% AUTHOR: h.schwab@tue.nl (Feb22)

function imgUS = getUSImage(receiveDataRFUS, fs,x_elem,c0, z,x)

    % manually set additional transducer settings:
    elemPos = [zeros(128,1), x_elem(:)];
    pulseDelays = zeros(128,1);
    fc = 8e6;
    
    % create Recon object and set paramters
    R = Recon2D(fs,elemPos,pulseDelays,c0,fc);
    R.reconParams.applyCompounding = false;
    R.visParams.dynamicRange       = 40;
    R.visParams.figHandle          = [];
    
    % initialize US image:
    imgUS = zeros(length(z),length(x),size(receiveDataRFUS,3));
    
    % loop over all frames:
    tic;
    for i = 1:size(receiveDataRFUS,3)
        
%         % retrieve reconstructed image:    
%         [RFData,dt]  = R.filter(receiveDataRFUS(:,:,i));
%         imgData      = R.reconMex(RFData, dt, z,x); 
%         imgUS(:,:,i) = R.visualize(imgData, z,x);
        
        % retrieve reconstructed image:         
        imgUS(:,:,i) = R.run(receiveDataRFUS(:,:,i),z,x);
        
        % print total computation time:
        toc;
    end

end