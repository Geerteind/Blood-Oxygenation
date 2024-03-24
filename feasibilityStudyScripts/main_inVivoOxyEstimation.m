%% MAIN_ININVIVOOXYESTIMATION oxygenation estimation of in vivo study 
%
% DESCRIPTION: 000 FILL IN A DOCUMENTATION 000
%
%AUTHOR: 000 TEAM 5 000

addpath('helpers');
s = settings;
        s.matlab.appearance.figure.GraphicsTheme.TemporaryValue= 'light'; % Can be set to auto (default), light, or dark
%% Settings

%  set file properties (note that W1->850nm and W2->750nm):
% Given data:
%filePathUS   = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\feasibilityStudyScripts\feasibilityStudyScripts\data\hans_arm_invivo_000_Rf_021722_174806_OBP_B_Extract_257.raw"; % string with path to raw file of ultrasound data to be loaded
%filePathW750 = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\feasibilityStudyScripts\feasibilityStudyScripts\data\hans_arm_invivo_000_Rf_021722_174806_OBP_PA_Wave1_PCAVG_384_128.raw"; % string with path to raw file of photoacoutic data (second wavelength) to be loaded
%filePathW850 = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\feasibilityStudyScripts\feasibilityStudyScripts\data\hans_arm_invivo_000_Rf_021722_174806_OBP_PA_Wave2_PCAVG_384_127.raw"; % string with path to raw file of photoacoutic data (first wavelength) to be loaded

% In vivo data: 
filePathUS   = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\InVivoData_group5\raw\PVA_blood_2_75MHz_min20deg_074_Rf_120523_211045_OBP_B_Extract_68.raw"; % string with path to raw file of ultrasound data to be loaded
filePathW750 = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\InVivoData_group5\raw\PVA_blood_2_75MHz_min20deg_074_Rf_120523_211045_OBP_PA_64_409.raw"; % string with path to raw file of photoacoutic data (second wavelength) to be loaded
filePathW850 = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\InVivoData_group5\raw\PVA_blood_2_75MHz_min20deg_074_Rf_120523_211045_OBP_PA_PCAVG_384_68.raw"; % string with path to raw file of photoacoutic data (first wavelength) to be loaded
%  set medium properties:
c0             = 1540; % scalar with speed of sound use din reconstruction [m/s]
mua_HBO2       = [2.7738,5.6654]; % Nwl element vector with the absorption coefficients of oxygenated blood @ [750,850]nm [a.u.]
mua_HB         = [7.5248,3.7019]; % Nwl element vector with the absorption coefficients of oxygenated blood @ [750,850]nm [a.u.]
mua_water      = [0.028484,0.041986]./100; % absorption coefficient of water @ [750,850]nm [1/m]     % plot(lambdaH2O,muH2O)
mua_background = [0.4,0.15]./100; % absorption coefficient of general tissue @ [750,850]nm [1/m] %      plot(lambda_genTiss,mu_abs_genTiss_30water)
mus_background = [11,9]./100;     % scattering coefficient @[750,850]nm of the background medium [1/m] 

% ... 000 add more properties if you need them ...

%  set image properties:
resolution_x_m  = 0.0001; % spacing between pixels in horizontal (x) direction [m]
resolution_z_m  = 0.0001; % spacing between pixels in vertical (z) direction [m]
imgBoundary_x_m = 0.0212; % distance of center of the transducer to left and right image boundary [m]
imgBoundary_z_m = 0.02; % distance of transducer surface to bottom image boundary [m]

%  set reconstruction properties:
framesIndexes = 30:44; % vector with indexes of frames that will be averaged (use the script "test_findFrames.m" to identify which frames to use!)
FNumber       = 2;   % positive number that defines the ratio of pixel depth to active aperture

%% Load receive data

%  load dataset from raw file (including metadata):
[receiveDataRFUS , fs_us,x_elem] = loadDataAcousticX(filePathUS,'US');
[receiveDataRF850, fs_pa       ] = loadDataAcousticX(filePathW850,'PA');
receiveDataRF750                 = loadDataAcousticX(filePathW750,'PA');

%% Precomputation

%  compute image axes according to settings: 
%  (with these two vectors, each image pixel gets assigned a location in 2D space)
z_axis    = (               0 : resolution_z_m : imgBoundary_z_m)'; % position of all pixels in z [m]
x_axis    = (-imgBoundary_x_m : resolution_x_m : imgBoundary_x_m)  ; % position of all pixels in x [m]

%  get lengths of axes:
Nx = length(z_axis);
Nz = length(z_axis);

%% Get ultrasound image:
%  In this section, the ultrasound image is retrieved. It serves as
%  geometric refecence, since morphological features can hardly be identified
%  in PA images. You can decide, which frame to retrieve and later plot an
%  overlay of PA and US (see the function "imgOverlay") if you like.

%  get ultrasound image:
figure(2); 
i_frameUS = 25; % index of frame to retrieve
imgUS = getUSImage(receiveDataRFUS(:,:,i_frameUS), fs_us,x_elem,c0, z_axis,x_axis);
% Determine water in US image
tissueMask = (imgUS > -35);

% Calculates minimum amount of water pixels for cropping
minWater = Nz; % Minimum amount of water pixels in a vertical line
for x = 1:(Nx-1)
    if (minWater > min(sum(~tissueMask(:,x)),sum(~tissueMask(:,x+1))))
        minWater = min(sum(~tissueMask(:,x)),sum(~tissueMask(:,x+1)));
    end
end
minWater = minWater - 5;
USsmall= imgUS(minWater:end,:); %The cropped US image
z_axis_cropped = z_axis(minWater:end);
%  show images:
subplot(121);   imagesc(x_axis*1e3,z_axis*1e3,imgUS); colormap gray; axis image; colorbar; % The normal US image
                %yline(z_axis(end)-minWater,'-','Cropping boundary');
                xlabel('x [mm]'); ylabel('z [mm]'); title('Image (US)');
subplot(122);   imagesc(x_axis(:)*1e3,z_axis_cropped*1e3,USsmall); colormap gray; axis image; colorbar; % The cropped US image
                xlabel('x [mm]'); ylabel('z [mm]'); title('Image (US) cropped');
                drawnow;

%% Apply frame averaging
%  Apply averaging of the receive data by computing the mean of multiple
%  frames (over the third dimension of hte array). Use only the frames defined 
%  in "framesIndexes" to be included in the averaging. Note that you can use 
%  the script "test_findFrames" to identify which frames should be averaged.

%  get mean of frames defined in "framesIndexes": 
% (convert Nz-by-Nx-by-Nfr array to a Nz-by-Nx array)
receiveDataRF750_avg = mean(receiveDataRF750(:,:,framesIndexes), 3);
receiveDataRF850_avg = mean(receiveDataRF850(:,:,framesIndexes), 3);

%% Reconstruct PA images 
%  reuse and adapt your code from the product development phase!

%  run PA reconstruction for both wavelengths (Nz-by-Nx-by-Nfr array):
imgDataRF750 = applyPAReconstruction(receiveDataRF750_avg,fs_pa,c0,x_elem,z_axis,x_axis,FNumber);
imgDataRF850 = applyPAReconstruction(receiveDataRF850_avg,fs_pa,c0,x_elem,z_axis,x_axis,FNumber);

% show PA images:
figure(1);
subplot(121); imagesc(x_axis*1e3, z_axis*1e3, imgDataRF750); 
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('Image @750nm');                                                        
subplot(122); imagesc(x_axis*1e3, z_axis*1e3, imgDataRF850);
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('Image @850nm');
drawnow;

%% Image postprocessing
%  Before you can apply the spectral unmixing, the image data requires some
%  further processing, which is done in this section. There might be more
%  processing steps needed, than for the simulation data, such as smoothing, 
%  thresholding, cropping or other post processing steps. Be creative! 

% Visualize the data distribution of the images before post-processing
figure(3);
subplot(121); histogram(imgDataRF750)
              axis image; 
              xlabel('pixel value'); ylabel('Frequency'); title('Data distribution image @750nm');
subplot(122); histogram(imgDataRF850)
              axis image; 
              xlabel('pixel value'); ylabel('Frequency'); title('Data distribution image @850nm');
drawnow;
%% Actual processing
%  apply envelope detection (Nz-by-Nx array):
f1 = 20;
imgDataComp750 = envelope(imgDataRF750,f1,'analytic'); % signal envelope using hilbert filter
imgDataComp850 = envelope(imgDataRF850,f1,'analytic');

%  correct for laser intensity ratio:
imgDataComp750 = imgDataComp750 * 1.63; %(laser intensity of 750 is 1.63 times lower than 850)

%  calculate fluence compensation weights (Nz-by-Nx array):
%  (you can segment the US image to assing materials, especially to 
%  distinguish tissue and the water layer on top of the tissue)
fluenceCompMap750 = calculateFluenceCompensationMap(z_axis,x_axis, [mua_water(1),mua_background(1)], imgUS);
fluenceCompMap850 = calculateFluenceCompensationMap(z_axis,x_axis, [mua_water(2),mua_background(2)], imgUS);

%  apply fluence compensation (Nz-by-Nx array):
imgDataComp750 = imgDataComp750 .* fluenceCompMap750; 
imgDataComp850 = imgDataComp850 .* fluenceCompMap850; 
%  cropping:
%(set regions that contain artifacts (e.g. at the top of the images) to 0)

imgDataComp750 = imgDataComp750(minWater:end,:);
imgDataComp850 = imgDataComp850(minWater:end,:);
%  smoothing:
imgDataComp750 = imgaussfilt(imgDataComp750,1) ; % Gaussian filter: maybe a filter that retain edges would be better?
imgDataComp850 = imgaussfilt(imgDataComp850,1) ; 

% Alternative smoothing: min filter then median filter
%se = strel('cube',2);
%imgDataComp750 = imerode(imgDataComp750,se);
%imgDataComp850 = imerode(imgDataComp850,se);
%imgDataComp750 = medfilt2(imgDataComp750,[5,5]);
%imgDataComp850 = medfilt2(imgDataComp850,[5,5]);

%  other processing steps:
% ... 000 ...
%%
% Visualize the data distribution of the images after post-processing
figure(3);
subplot(121); histogram(imgDataComp750)
              xlabel('pixel value'); ylabel('Frequency'); title('Data distribution image @750nm');
subplot(122); histogram(imgDataComp850)              
              xlabel('pixel value'); ylabel('Frequency'); title('Data distribution image @850nm');
drawnow;
%% Final images after processing
%  show images:
figure(3);
subplot(121); imagesc(x_axis*1e3, z_axis_cropped*1e3, imgDataComp750); 
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('Processed image @750nm');
subplot(122); imagesc(x_axis*1e3, z_axis_cropped*1e3, imgDataComp850); 
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('Processed image @850nm');
drawnow;

%% Concentration estimation
%  reuse and adapt your code from the product development phase!

[imgDataHB,imgDataHBO2] = applySpectalUnmixing(imgDataComp750,imgDataComp850, mua_HB,mua_HBO2);

% show images:
figure(5); 
subplot(121); imagesc(x_axis*1e3, z_axis_cropped*1e3, imgDataHB); 
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('c_{HB}');
subplot(122); imagesc(x_axis*1e3, z_axis_cropped*1e3, imgDataHBO2); 
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('c_{HBO2}');
drawnow;
%% Oxygenation estimation              
%  reuse and adapt your code from the product development phase!
%  convert component data (Nz-by-Nx array) into oxygenation map (Nz-by-Nx array) [%]:
oxyMap = 100*(imgDataHBO2./(imgDataHB+imgDataHBO2));
% define a boolean matrix (Nz-by-Nx) that is true for each pixel, for 
%   which the image value is lower than a 40% of the maximum in the unmixed
%   images
zeroMask = ( imgDataComp750+imgDataComp850 < 0.3*max(max(imgDataComp750(:)),max(imgDataComp850(:))) )...
          |( oxyMap > 100 );
oxyMap(zeroMask) = 0;
oxyMap(oxyMap<0) = 0;
%  show image:
lowerThresh = 0.05; % lowest oxygenation value that will be displayed on colormap
figure(6); imagesc(x_axis*1e3, z_axis_cropped*1e3, oxyMap, [lowerThresh,100]); 
           colormap([zeros(1,3); cool]); axis image; colorbar; 
           xlabel('x [mm]'); ylabel('z [mm]'); title('oxygenation [%]');
%%
overlayOxy = imgOverlay(USsmall,oxyMap,0);
figure(7); imagesc(x_axis*1e3, z_axis_cropped*1e3, overlayOxy, [lowerThresh,100]); 
           colormap([zeros(1,3); hot]); axis image; colorbar; 
           xlabel('x [mm]'); ylabel('z [mm]'); title('oxygenation [%] with US image overlayed for structural information');

