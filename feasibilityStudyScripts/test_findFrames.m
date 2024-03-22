%% TEST_FINDFRAMES reconstructs all PA frames and helps finding out which frames to average
% in the main script, the receive data are averaged, such that only one
% image per wavelengths needs to be reconstructed. Since the reconstruction
% is a linear process, it is mathematically the same to reconstruct all
% images first and to average them after. This way you can flexibly find
% out, which frames should be used for the averaging in the main script.
% The main aspect of this search is to find a good trade-off between noise
% decorrelation and signal decorrelation due to motion.
% The findings can be used to set the parameter "framesIndexes" in the 
% main script after!
%
%AUTHOR: h.schwab@tue.nl (Feb22)
addpath('helpers');

%% Settings:

% set paths (note that W1->850nm and W2->750nm):
%absoluteFilePath = 'C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\feasibilityStudyScripts\feasibilityStudyScripts\';
%filePathUS   = append(absoluteFilePath,'data/hans_arm_invivo_000_Rf_021722_174806_OBP_B_Extract_257.raw'); % string with path to raw file of ultrasound data to be loaded
%filePathW750 = append(absoluteFilePath,'data/hans_arm_invivo_000_Rf_021722_174806_OBP_PA_Wave2_PCAVG_384_127.raw'); % string with path to raw file of photoacoutic data (second wavelength) to be loaded
%filePathW850 = append(absoluteFilePath,'data/hans_arm_invivo_000_Rf_021722_174806_OBP_PA_Wave1_PCAVG_384_128.raw'); % string with path to raw file of photoacoutic data (first wavelength) to be loaded
% In vivo data: 
filePathUS   = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\InVivoData_group5\raw\PVA_blood_2_75MHz_min20deg_078_Rf_120523_213855_OBP_B_Extract_36.raw"; % string with path to raw file of ultrasound data to be loaded
filePathW750 = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\InVivoData_group5\raw\PVA_blood_2_75MHz_min20deg_078_Rf_120523_213855_OBP_PA_64_221.raw"; % string with path to raw file of photoacoutic data (second wavelength) to be loaded
filePathW850 = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\InVivoData_group5\raw\PVA_blood_2_75MHz_min20deg_078_Rf_120523_213855_OBP_PA_PCAVG_384_36.raw"; % string with path to raw file of photoacoutic data (first wavelength) to be loaded
% set image properties:
resolution_x_m  = .2*1e-3; % spacing between pixels in horizontal (x) direction [m]
resolution_z_m  = .1*1e-3; % spacing between pixels in vertical (z) direction [m]
imgBoundary_x_m = 20*1e-3; % distance of center of the transducer to left and right image boundary [m]
imgBoundary_z_m = 15*1e-3; % distance of transducer surface to bottom image boundary [m]

% set reconstruction properties:
c0      = 1500; % scalar speed of sound for reconstruction [m/s]
FNumber = 2;   % positive number that defines the ratio of pixel depth to active aperture

%% Load all frames:

% load dataset from raw file (including metadata):
[receiveDataRFUS , fs_us,x_elem] = loadDataAcousticX(filePathUS,'US');
[receiveDataRF850, fs_pa       ] = loadDataAcousticX(filePathW850,'PA');
receiveDataRF750                 = loadDataAcousticX(filePathW750,'PA');

% precompute image axes:
z_axis    = (               0 : resolution_z_m : imgBoundary_z_m)'; % position of all pixels in z [m]
x_axis    = (-imgBoundary_x_m : resolution_x_m : imgBoundary_x_m)  ; % position of all pixels in x [m]

%% Reconstruct all frames:

% run US reconstruction:
imgsDataUS = getUSImage(receiveDataRFUS, fs_us,x_elem,c0, z_axis,x_axis);

% run PA reconstruction for both wavelengths (Nz-by-Nx-by-Nfr array):
imgsDataRF750 = applyPAReconstruction(receiveDataRF750, fs_pa,c0, x_elem, z_axis,x_axis,FNumber);
imgsDataRF850 = applyPAReconstruction(receiveDataRF850, fs_pa,c0, x_elem, z_axis,x_axis,FNumber);

% save data:
%(you can uncommend this part to avoid the need to reconstruct them again)
save('ReconstructedRDData.mat','z_axis','x_axis','imgsDataRF750','imgsDataRF750','imgsDataUS');

%% View all frames:

% set indexes of frames that will be displayed:
%(adjust this value to find a good time span where the image content does not move to much)
%framesIndexes =  1:size(imgsDataRF750,3);% 5:15;%  vector with indexes of frames that will be averaged 
framesIndexes =  15:27; % 5:15;%  vector with indexes of frames that will be averaged 
for i_fr = framesIndexes
    
    % show US image:
    figure(1); imagesc(x_axis*1e3,z_axis*1e3,imgsDataUS(:,:,i_fr));
               colormap gray; axis image; colorbar; 
               xlabel('x [mm]'); ylabel('z [mm]'); title(['Image ',num2str(i_fr),' (US)']);
    
    % show PA images:
    figure(2);
    subplot(121); imagesc(x_axis*1e3,z_axis*1e3,abs(hilbert(imgsDataRF750(:,:,i_fr)))); 
    % subplot(121); imagesc(x_axis*1e3,z_axis*1e3,imgOverlay(imgUS(:,:,i_fr),abs(hilbert(imgDataRF750(:,:,i_fr))),.05,.2)); 
                  colormap hot; axis image; colorbar; 
                  xlabel('x [mm]'); ylabel('z [mm]'); title(['Image ',num2str(i_fr),' @750nm']);                                                        
    subplot(122); imagesc(x_axis*1e3,z_axis*1e3,abs(hilbert(imgsDataRF850(:,:,i_fr))));
    % subplot(122); imagesc(x_axis*1e3,z_axis*1e3,imgOverlay(imgUS(:,:,i_fr),abs(hilbert(imgDataRF850(:,:,i_fr))),.05,.2));
                  colormap hot; axis image; colorbar; 
                  xlabel('x [mm]'); ylabel('z [mm]'); title(['Image ',num2str(i_fr),' @850nm']);
                     
   % pause (set display rate by argument [s]):
   pause(.01);
    
end

%% Apply averaging:

framesIndexes =  15:27;

% apply frame averaging on reconstructed images:
imgDataRF750_avg = mean(imgsDataRF750(:,:,framesIndexes), 3);
imgDataRF850_avg = mean(imgsDataRF850(:,:,framesIndexes), 3);

% show averaged images:
figure(2);
subplot(121); imagesc(x_axis,z_axis,abs(hilbert(imgDataRF750_avg))); 
              colormap hot; axis image; title('PA image @750nm');
subplot(122); imagesc(x_axis,z_axis,abs(hilbert(imgDataRF850_avg))); 
              colormap hot; axis image; title('PA image @850nm');
 

