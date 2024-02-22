% Analysis4_atHome
% This script shows the bandlimiting bahavior of ultrasound, how it affects
% the image appearance and what can be done to display an image that still
% looks good.
%
% Instruction: 
% Read the caption of each section and then run the respective section
% using CTRL+ENTER.
%
% NOTE: please add the subfolder "helpers" and its subfolders to your
% matlabpath!
%

%% create vessel phantom:
% If we had constant fluence in the tissue, the cross section of a vessel
% would look like a disk:

% create image axes [m]:
z = (      0: 0.1e-3 : 40e-3)';
x = (-20e-3 : 0.1e-3 : 20e-3);

% create phantom as disk at center of the image with 5mm radius:
p0 = zeros(length(z),length(x));
p0(sqrt((z-20e-3).^2 + x.^2) < 5e-3) = 1;

% plot phantom:
figure(1); subplot(221); imagesc(z*1e3,x*1e3,p0); title('Initial pressure image');
axis image; colormap gray; colorbar; xlabel('x [mm]'); ylabel('z [mm]');

%% apply point-spread function:
% A simple model to describe the effect of the imaging system on the image
% content in photoacoustic imaging (measurement plus reconstruction) is to
% convolve the original phantom with a point spread function (PSF). Here, 
% an approximate PSF is computed as a 2D kernel that consists of a
% windowed sine function in vertical (axial) direction multiplied with a 
% Gaussian function in horizontal (lateral) direction.
% How would oyu describe the resulting radio-frequency (RF) image compared 
% to theinitial initial pressure distribution?

% create simplified point-spread function:
psf = .15 * gauspuls(-10:10,.3)' .* exp(-.01*(-20:20).^2);

% plot filter kernel:
figure(1); subplot(222); imagesc((-20:20)*.1,(-10:10)*.1,psf); title('point-spread function');
axis image; colormap parula; colorbar; xlabel('x [mm]'); ylabel('z [mm]');

% filter phantom with psf to get a simple photacoustic raw-data (RF) image:
imgRF = imfilter(p0,psf);

% plot RF image:
figure(1); subplot(223); imagesc(z*1e3,x*1e3,imgRF); title('RF image');
axis image; colormap gray; colorbar; xlabel('x [mm]'); ylabel('z [mm]');

%% envelope detection:
% The RF image we have now is the unprocessed image, a photoacoutic machine
% will create. The strong oscillations in axial direction are caused by the 
% limited frequency range (bandwidth) of the transducer. This leads to an 
% edge filtering effect, such that we only asee a signal on the edges of the 
% disk. This osciallating image content is difficult to interpret with our
% eyes. For that reason, we usually apply an envelope detection 
% (demodulation), which will essentially connect the peaks of the signal to
% get a smooth image with strictly positive values.
% Therefore, the image you can expect of a vessel rather resembles two dots
% above each other than a full circle. Still, the intensity of these edges 
% can be used to estimate the oxygenation within the vessel!   

% apply envelope detection:
imgEnvelope = abs(hilbert(imgRF));

% plot envelope image:
figure(1); subplot(224); imagesc(z*1e3,x*1e3,imgEnvelope); title('Envelope image');
axis image; colormap gray; colorbar; xlabel('x [mm]'); ylabel('z [mm]');

% plot horizontal and lateral line profiles:
figure(2); 
subplot(211); plot(z*1e3,p0(:,200), z*1e3,imgRF(:,200), z*1e3,imgEnvelope(:,200))
    title('axial profile'); xlabel('z [mm]'); legend('p0','RF','Envelope');
subplot(212); plot(x*1e3,p0(200,:), x*1e3,imgRF(200,:), x*1e3,imgEnvelope(200,:)); 
    title('lateral profile'); xlabel('z [mm]'); legend('p0','RF','Envelope')


