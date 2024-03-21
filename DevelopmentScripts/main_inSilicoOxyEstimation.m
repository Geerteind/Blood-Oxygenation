
%% MAIN_INSILICOOXYESTIMATION template for oxygenation estimation in simulation data 
%
%DESCRIPTION:
% This is the main script for your software development work. Please
% download the simulation dataset for your team from canvas, put it into 
% the "data" folder and start to code! 
% In this script, we propose a general structure for your code that you can 
% build upon by filling in the missing parts (000). But you can also 
% freely change everything to your own liking and include any additional 
% processing steps you can think of!
%
%SOME TIPS:
%  - run: ">>doc FUNCTIONNAME" if you don't know how to use a function! 
%  - use CTR+ENTER to run a single section and F9 to run only the code you highlighted
%  - if any function you write needs more information you can easily define 
%    more input or output arguments!
%  - always take care of using consistent units, generally work with [m],[s]
%  - if you cannot find out what goes wrong, plot all the quantities you 
%    computed using commands, such as: "imagesc", "plot", "scatter" 
%
%GEOMETRY:
% This is the geometry convention used thoughout the script:
%
%    o- - - - -> x[m]
%    |\ 
%    | \ theta [deg] 
%    |_/\   (-180,180)
%    |   \
%    v
%   z[m]
%
%  (in all arrays, the z position is the first entry, the x position is the second)
%
% HAVE FUN !
%
%AUTHOR: 000 TEAM 5 000

addpath('helpers');
%% Settings
%  in this section, you can define any parameters that are used throughout
%  the script. If you want to add any tunable value, you can add it here to
%  keep a nice overview!

%  set path to matfile:
filePath = "C:\Users\annel\MATLAB Drive\Blood-Oxygenation\DevelopmentScripts\data\Data_group5.mat" ; % string with path to raw file of ultrasound data to be loaded
%filePath = "C:\Users\annel\OneDrive - TU Eindhoven\Year 3\Q3\DBL Blood oxygenation\PhantomExperiment\raw\PVA_blood_2_75MHz_min20deg_071_Rf_112323_231727_OBP_B_1536.raw"

%  set medium properties 
%  (according to the values that were posted in the anonuncement on canvas):
mua_HBO2       = [2.7738,5.6654] ; % 2 element vector with the absorption coefficients @(750,850)nm of oxygenated blood [1/cm]
mua_HB         = [7.5248,3.7019] ; % 2 element vector with the absorption coefficients @(750,850)nm of deoxygenated blood [1/cm]
g_blood        = 0.95 ; % scalar anisotropy factor of both kinds of blood [1]
mua_background = [0.4,0.15]./100 ; % absorption coefficient @[750,850]nm of the background medium [1/cm]
mus_background = [11,9]./100 ; % scattering coefficient @[750,850]nm of the background medium [1/cm] 
g_background   = 0 ;  % scalar anisotropy factor of background medium [1]

%  set image properties:
resolution_x_m  = 0.0001 ; % spacing between pixels in horizontal (x) direction [m]
resolution_z_m  = 0.0001 ; % spacing between pixels in vertical (z) direction [m]
imgBoundary_x_m = 0.025 ; % distance of center of the transducer to left and right image boundary [m]
imgBoundary_z_m = 0.05 ; % distance of transducer surface to bottom image boundary [m]

%  set reconstruction properties:
FNumber = 2; % positive number that defines the ratio of pixel depth to active aperture

%%  Load receive data
% Load the simulated receive data sets. Run:
%   >> doc loadDataSimulation
% to understand the inputs and outputs of the function, if required.

%  load dataset from mat file (including metadata):
[receiveDataRF750, receiveDataRF850, fs, c0, x_elem] = loadDataSimulation(filePath);

%% Precomputation:
%  Based on the settings you defined before, some variables need to be
%  precomputed before starting the actual data processing. .

%  calculate effective light attenuation coefficients in the backgroud medium:
mueff_background = sqrt(3*mua_background.*(mua_background+(1-g_background)*mus_background)) ; 

%  compute image axes according to settings: 
%  (with these two vectors, each image pixel gets assigned a location in 2D space)
z_axis    = (               0 : resolution_z_m : imgBoundary_z_m)'; % position of all pixels in z [m]
x_axis    = (-imgBoundary_x_m : resolution_x_m : imgBoundary_x_m)  ; % position of all pixels in x [m]

%% Reconstruct PA images 
% In this section, the receive data is reconstructed to get actual image data.
% For this, you need to fill in the essential parts missing in the function
% "applyPAReconstruction" using the information from above. You can tune
% the image quality / resolution by changing the value for the FNumber,
% which defines how large the active receive aperture is. Run:
%   >> edit applyPAReconstruction
% to start writing your code according to what you have found out in the 
% Research phase! 

%  run PA reconstruction for both wavelengths (Nz-by-Nx array):
imgDataRF750 = applyPAReconstruction(receiveDataRF750, fs,c0,x_elem, z_axis,x_axis,FNumber);
imgDataRF850 = applyPAReconstruction(receiveDataRF850, fs,c0,x_elem, z_axis,x_axis,FNumber);

%  show PA images:

figure(1);
subplot(121); imagesc(x_axis*1e3, z_axis*1e3, imgDataRF750); 
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('Image @750nm');                                                        
subplot(122); imagesc(x_axis*1e3, z_axis*1e3, imgDataRF850);
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('Image @850nm');
drawnow;

%%  Image postprocessing
%  Before you can apply the spectral unmixing, the image data requires some
%  further processing, which is done in this section. 

%  apply envelope detection (Nz-by-Nx array):
imgDataComp750 = logcomp(imgDataRF750) ; 
imgDataComp850 = logcomp(imgDataRF850) ; 

%  calculate fluence compensation weights (Nz-by-Nx array):
%  (calculate a depth dependent correction map according to what you found 
%   out in the Research phase. Remeber that the axis is in [m] and mu is in [1/cm])
fluence_mua = mua_background; 
fluence_mus = mus_background;

fluenceCompMap750 = 1./exp(-z_axis*mueff_background(1)) ; % used the fluence formula
fluenceCompMap850 = 1./exp(-z_axis*mueff_background(2)) ; 

%  apply fluence compensation (Nz-by-Nx array):
imgDataComp750 = fluenceCompMap750 .* imgDataComp750 ; 
imgDataComp850 = fluenceCompMap850 .* imgDataComp850 ; 

%  apply smoothing if you want:
%  (search the internet of how you can smooth an image using filter functions in Matlab)
%imgDataComp750 = 000 ; 
%imgDataComp850 = 000 ; 

%  show images:
figure(3);
subplot(121); imagesc(x_axis*1e3, z_axis*1e3, imgDataComp750); 
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('Processed image @750nm');
subplot(122); imagesc(x_axis*1e3, z_axis*1e3, imgDataComp850); 
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('Processed image @850nm');
drawnow;

%  show fluence compensation maps:
figure(4);
maxVal = max(max(fluenceCompMap750(:)),max(fluenceCompMap850(:))); %(maximum of both  to scale plot)
subplot(121); imagesc(x_axis*1e3, z_axis*1e3, fluenceCompMap750, [1,maxVal]); 
              colormap jet; axis image; colorbar;
              xlabel('x [mm]'); ylabel('z [mm]'); title('Fluence compensation map @750nm');                                                       
subplot(122); imagesc(x_axis*1e3, z_axis*1e3, fluenceCompMap850, [1,maxVal]); 
              colormap jet; axis image; colorbar;
              xlabel('x [mm]'); ylabel('z [mm]'); title('Fluence compensation map @850nm');
drawnow;

%% Concentration estimation
% This section contains the actual spectral unmixing part. Here, the
% images at different wavelength (imgDataComp750 and imgDataComp850) are
% converted into images of different components (imgDataHB and imgDataHBO2).
% The unmixing is carried out in the function "applySpectalUnmixing", which 
% you need to fill with content. Run:
%   >> edit applySpectalUnmixing
% to start writing your code for that!
% Note that the concentrations in the component data can be in arbitrary
% units, because they will be set into relation afterwards and any scaling 
% will cancel out!

%  unmix wavelength data (Nz-by-Nx array) to component data (Nz-by-Nx array):
[imgDataHB,imgDataHBO2] = applySpectalUnmixing(imgDataComp750,imgDataComp850, mua_HB,mua_HBO2);
                                                                                                                                                       
%  show images:
figure(5); 
subplot(121); imagesc(x_axis*1e3, z_axis*1e3, imgDataHB); 
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('c_{HB}');
subplot(122); imagesc(x_axis*1e3, z_axis*1e3, imgDataHBO2); 
              colormap hot; axis image; colorbar; 
              xlabel('x [mm]'); ylabel('z [mm]'); title('c_{HBO2}');
drawnow;

%% Oxygenation estimation              
%  You now have maps of concentrations of the two components of blood and
%  only need to convert them (pixel-by-pixel) to an oxygenation map.
%  After this, it is advised to mask out (set to 0) areas, that do not
%  contain vessels, because an oxygenation value will be computed in each 
%  pixel, but the value only makes sense if the actual components in that 
%  pixel were oxygenated and deoxygenated blood! For this masking, you can 
%  set the oxygenation value to 0 for all pixels that did not have high 
%  image intensities in the original PA images.
%  To get your final oxygenation score at the location of the vessels, you 
%  can plot a lineplot of the oxygenation map ("plot") or you can use the 
%  "Data Tips" tool in the image Figure!

%  convert component data (Nz-by-Nx array) into oxygenation map (Nz-by-Nx array) [%]:
oxyMap = normalize(imgDataHBO2./(imgDataHB+imgDataHBO2)); 

%  define mask to 0 values:
%  (define a boolean matrix (Nz-by-Nx) that is true for each pixel, for 
%   which the image value is lower than a certain threshold):
zeroMask = ( oxyMap<0.1 )... 
          |( oxyMap == 0 );
%  apply that mask:
%  (by settingthe oxygenation map to 0 at all true pixels of the mask to 0)      
oxyMap(zeroMask) = 0;

%  show image:
lowerThresh = 0.05; % lowest oxygenation value that will be displayed on colormap
figure(6); imagesc(x_axis*1e3, z_axis*1e3, oxyMap, [lowerThresh,100]); 
           colormap([zeros(1,3); cool]); axis image; colorbar; 
           xlabel('x [mm]'); ylabel('z [mm]'); title('oxygenation [%]');

