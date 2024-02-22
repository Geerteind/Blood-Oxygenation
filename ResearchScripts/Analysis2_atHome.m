% Analysis2_atHome This script relies on aquired AcousticX PA data of a 
% single thread that moves over frames in time to find out the
% photoacoustic reconstruction delay! 
% 
% Instruction:
% - Add the subfolder "helpers" and its subfolders to your matlabpath!
% - Run the script using the run botton or press "F5".
% - You will be asked to chose a file. Please select the ".raw"-file you
%    saved during the experiments that contains photoacoustic receive data 
%    of a single thread  
% - Familiarize yourself with the two images on  the left (receive data)
%    and on the right(image data). Understand the axes and what they show!
% - Look at the image data (right side) and try to come up with an equation
%    that describes the travel time from the target to each trnasducer
%    element. 
% - Type the correct equation into the text field below the images as a
%    function of "x_el","x","z" and "c0" and press enter.
% - Verify your delay equation by comparing the resulting curve to the
%    intensity distribution in the receive data!
% - Close the window to end the gui!
%
% h.schwab@tue.nl (Feb23)

%% Settings:

% set the speed of sound that you found in the experiment:
c0 = 1480; % speed of sound [m/s]

%% Load AcousticX raw dataset:

% select raw file path:
[fileName,pathName] = uigetfile('*.raw','Select raw file containing PA data');
filePath = [pathName,filesep,fileName];
% or set manually:
% filePath = 'Arm_Min_006_Rf_021122_201111-03_OBP_PA_64_820.raw';

% load PA data:
[RXData, fs,x_elem] = loadDataAcousticX(filePath);
% delete laser distortion:
RXData(1:100,:,:) = 0;

%% generate RX Data:

% acquisition settings:
fc = fs/4; % center frequency [Hz]
da = abs(diff(x_elem(1:2)));
x = x_elem(1):da/2:x_elem(end); % grid vector in x-direction [m] 
z = 0:da/4:size(RXData,1)/fs*c0; % grid vector in x-direction [m] 
PALaserDelay = -15/fc; % PA laser pulse delay [s]
normalAnglesRX = zeros(length(x_elem),1); % normal angles of all elements relative to z-axis [deg]
RXPos = [zeros(length(x_elem),1), x_elem(:)];

% create reconstructor object:
R = Recon2D(fs,RXPos,PALaserDelay,c0,fc);
% R.setTXDelaysPlaneWave(TXAngles);

%% Window for DAS Recon (1):

% initialize figure
fig1 = figure(1); clf
fig1.Units = 'normalized';
fig1.Position = [.2,.1,.6,.8];

% adjust recon params:
R.visParams.supressDrawnow = true;
R.reconParams.apodCutoffAngle = 20;
R.visParams.figHandle = fig1;
R.visParams.dynamicRange = 100;

% define image axes:
t = (0:size(RXData,1)-1)/R.fs;
x_el  = R.elemPos(:,2)'*1e3;
nSkip = 4;

% create GUI elements:
ax1 = axes('Parent',fig1,'position',[0.1 0.3  0.35 0.6]);
ax2 = axes('Parent',fig1,'position',[0.6 0.3  0.35 0.6]);
b = uicontrol( 'Parent',fig1,'Style','edit','Units','normalized','Position',[.25,.1,.5,.05]...
              ,'String','1/z_pix.^2 + (x_el-x_pix).^2/c0', 'FontSize',12, 'FontName','Monospaced'); 
b2 = annotation( gcf,'textbox','String','change delay:   \tau =  ', 'FontSize',12,'units','normalized', 'Position', [.0,.1,.25,.05]...
                ,'Interpreter','Tex','HorizontalAlignment','right','EdgeColor',fig1.Color);
         
% mimick external function call in acquisition loop:
i=0;
while ishandle(fig1)
    
    % increase frame counter:
    i_fr = mod(i,size(RXData,3))+1;
    i=i+1;
            
    % update curve function from GUI text edit field to display in RXData:
    try
        curve = eval(b.String);
        curveCol = 'b'; 
    catch
        curve = ones(1,length(x_el));
        curveCol = 'r'; 
    end
%     disp(b.String);
    
    % reconstruct image:
%     imgDAS = R.reconPWCuda(RXData(:,:,i_fr), 1/R.fs, z, x, struct('apodAngle',20,'isPA',true));  
    imgDAS = R.reconMex(RXData(:,:,i_fr), 1/R.fs, z, x);  
    
    % find peak in image and set pixel pos to peak pos:
    imgSrchMax = abs(imgDAS);
    imgSrchMax(z*1e3<3,:,:) = 0;
    [z_pix, x_pix] = find(imgSrchMax==max(imgSrchMax(:)),1);
    x_pix = (x(1) + x_pix*(x(2)-x(1)))*1e3;
    z_pix = (z(1) + z_pix*(z(2)-z(1)))*1e3;
    
    % plot RXData image with curve plot:
    axes(ax1)
    imagesc(x*1e3,t,logcomp(RXData(:,:,i_fr),30)); axis normal; colormap gray; 
%     imagesc(x*1e3,ct*1e3,real(RXData(:,:,round(size(RXData,3)/2)))); axis image; colormap gray;
    ylim([0,t(end)]); title('Receive Data');
    ylabel('t [s]'); xlabel('x_el [mm]');
    hold on; plot(x_el,curve/1e3+PALaserDelay,curveCol,'LineWidth',2); hold off;
    hold on; plot([x_el(1:nSkip:end);x_el(1:nSkip:end)],[zeros(1,length(x_el(1:nSkip:end)));curve(1:nSkip:end)/1e3+PALaserDelay],'b','LineWidth',.1); hold off;
    
    % plot image data with pixel/element line plots:
    axes(ax2)
    imagesc(x*1e3,z*1e3,logcomp(imgDAS,30)); axis image; colormap gray;   
    ylim([0,z(end)]*1e3); title('Image Data');
    R.plotProbe; legend off;   
    ylabel('z [mm]'); xlabel('x [mm]');
    hold on; scatter(x_pix,z_pix, [curveCol,'o'], 'LineWidth',2); hold off;
    hold on; plot([x_el(1:nSkip:end);x_pix*ones(length(x_el(1:nSkip:end)))],[zeros(1,length(x_el(1:nSkip:end)));z_pix*ones(length(x_el(1:nSkip:end)))],'b','LineWidth',.1); hold off;
    text(x_pix,z_pix,' [x_pix,z_pix]','Color','w','FontName','Monospaced','FontSize',10,'interpreter','none','FontWeight','bold')
    
%     pause(1);
    drawnow;
    
end
