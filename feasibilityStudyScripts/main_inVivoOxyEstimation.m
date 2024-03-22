%% MAIN_ININVIVOOXYESTIMATION oxygenation estimation of in vivo study 
%
% DESCRIPTION: 000 FILL IN A DOCUMENTATION 000
%
%AUTHOR: 000 FILL IN TEAMNAME 000

addpath('helpers');

%% Settings

%  set file properties (note that W1->850nm and W2->750nm):
filePathUS   = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\feasibilityStudyScripts\feasibilityStudyScripts\data\hans_arm_invivo_000_Rf_021722_174806_OBP_B_Extract_257.raw"; % string with path to raw file of ultrasound data to be loaded
filePathW750 = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\feasibilityStudyScripts\feasibilityStudyScripts\data\hans_arm_invivo_000_Rf_021722_174806_OBP_PA_Wave1_PCAVG_384_128.raw"; % string with path to raw file of photoacoutic data (second wavelength) to be loaded
filePathW850 = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\Code\feasibilityStudyScripts\feasibilityStudyScripts\data\hans_arm_invivo_000_Rf_021722_174806_OBP_PA_Wave2_PCAVG_384_127.raw"; % string with path to raw file of photoacoutic data (first wavelength) to be loaded

%  set medium properties:
c0             = 1540; % scalar with speed of sound use din reconstruction [m/s]
mua_HBO2       = [2.7738,5.6654]; % Nwl element vector with the absorption coefficients of oxygenated blood @ [750,850]nm [a.u.]
mua_HB         = [7.5248,3.7019]; % Nwl element vector with the absorption coefficients of oxygenated blood @ [750,850]nm [a.u.]
mua_water      = 000; % absorption coefficient of water @ [750,850]nm [1/cm]     % plot(lambdaH2O,muH2O)
mua_background = [0.4,0.15]./100; % absorption coefficient of general tissue @ [750,850]nm [1/m] %      plot(lambda_genTiss,mu_abs_genTiss_30water)
mus_background = [11,9]./100;     % scattering coefficient @[750,850]nm of the background medium [1/m] 

% ... 000 add more properties if you need them ...

%  set image properties:
resolution_x_m  = 0.0001; % spacing between pixels in horizontal (x) direction [m]
resolution_z_m  = 0.0001; % spacing between pixels in vertical (z) direction [m]
imgBoundary_x_m = 0.0212; % distance of center of the transducer to left and right image boundary [m]
imgBoundary_z_m = 0.02; % distance of transducer surface to bottom image boundary [m]

%  set reconstruction properties:
framesIndexes = 000; % vector with indexes of frames that will be averaged (use the script "test_findFrames.m" to identify which frames to use!)
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
i_frameUS = 10; % index of frame to retrieve
imgUS = getUSImage(receiveDataRFUS(:,:,i_frameUS), fs_us,x_elem,c0, z_axis,x_axis);

%  show images:
imagesc(x_axis*1e3,z_axis*1e3,imgUS); colormap gray; axis image; colorbar;
xlabel('x [mm]'); ylabel('z [mm]'); title('Image (US)');
drawnow;

%% Apply frame averaging
%  Apply averaging of the receive data by computing the mean of multiple
%  frames (over the third dimension of hte array). Use only the frames defined 
%  in "framesIndexes" to be included in the averaging. Note that you can use 
%  the script "test_findFrames" to identify which frames should be averaged.

%  get mean of frames defined in "framesIndexes": 
% (convert Nz-by-Nx-by-Nfr array to a Nz-by-Nx array)
receiveDataRF750_avg = 000;
receiveDataRF850_avg = 000;

%% Reconstruct PA images 
%  reuse and adapt your code from the product development phase!

%  run PA reconstruction for both wavelengths (Nz-by-Nx-by-Nfr array):
imgDataRF750 = 000;
imgDataRF850 = 000;

% show PA images:
% ... 000 ...

%% Image postprocessing
%  Before you can apply the spectral unmixing, the image data requires some
%  further processing, which is done in this section. There might be more
%  processing steps needed, than for the simulation data, such as smoothing, 
%  thresholding, cropping or other post processing steps. Be creative! 

%  apply envelope detection (Nz-by-Nx array):
imgDataComp750 = 000;
imgDataComp850 = 000;

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
% ... 000 ...

%  smoothing:
% ... 000 ...

%  other processing steps:
% ... 000 ...

%  show images:
% ... 000 ...

%% Concentration estimation
%  reuse and adapt your code from the product development phase!

% ... 000 ...

%% Oxygenation estimation              
%  reuse and adapt your code from the product development phase!

% ... 000 ...


