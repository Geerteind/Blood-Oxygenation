% RECON2D Reconstruction class for 2D ultrasound imaging for various imaging 
% configurations
%
%NOTE: 
% This class was renamed from "DASReconstructor" in June 2020
%
%DESCRIPTION:
% This ultrasound reconstruction class aims to be as general as possible.
% It can be applied to any configuration of transducer elements and
% arbitrary wavefronts, defined only by the transmission delays of the
% elements (which can also be derived from other parameters using converter
% method). The class allows for GPU processing by using gpu arrays as 
% input. The core of the reconstruction is written as a function that 
% greatly benefits from parallelization. The class also comes with various
% data format and coordinate system converters, visualization tools to show
% the geometry, the filter spectra and the final US image, as well as a
% simple simulation tool to generate measurement data based on a set of
% arbitrary point scatterers. There are three data processing methods, one
% for filtering the receive data, which can also convert into IQ data if
% demanded, one for the actual reconstruction, where multiple options for
% space-dependent apodization and insonification weighting can be chosen
% from, and finally a visulaization method to display the reconstructed
% dataset. Furthermore, The method .run() that combines all three
% processing steps. The parameters for all three processing steps can be 
% set manually or by using a profile. 
%
%GEOMETRY:
%
%    o- - - - -> x[m]
%    |\ 
%    | \ theta [deg] 
%    |_/\   (-180,180)
%    |   \
%    v
%   z[m]
%
% (in all position descriptions, the z position is the first entry)
%
%EXAMPLES:
%
% create reconstructor object:
% >> dasR = Recon2D(fs, elemPos, TXPulseDelays, c0, fc);
% change a general settings:
% >> dasR.printCompTimes = true;
% change a filter parameter:
% >> dasR.filtParams.filtType = 'simpleBP';
% change a reconsruction parameter:
% >> dasR.reconParams.apodCutoffAngle = 15;
% change a visualization parameter:
% >> dasR.visParams.dynamicRange = 70;
% run all three processing steps (filter, reconstruct, visualize):
% >> dasR.run(RXData,z,x);
% (Please also refer to the example macro to learn more about the 
%  functionality of this class)
%
%PROPERTIES:
%
%Required Recon Properties:
% .fs               sampling frequency [1/s]
% .elemPos          Nx2 array with centerPositions of the elements [m] in order [zPos,xPos]
% .TXPulseDelays    (NaTX x Nacq)-matrix with delays of the TX pulse peak for each element [s]
% .c0               scalar speed of sound [m/s]
%
%Optional Recon-Related Properties:
% .fc               center frequency [1/s] (not mandatory but required for some filters)
% .normalAngleRX    .Na-element vector with normal angles of the elements 
%                   (or scalar if equal for all elements). Default: 0
% .elemPosTX        Nx2 array with centerPositions of the TX elements (if different) [m] in order [zPos,xPos] 
%
%Other optional properties:
% .elemWidth       scalar containing the width of each transducer element [m]
% .bandWidth       scalar containing fractional bandwidth relative to .fc [1] (decimal, not in % !)
%
%Parameter Settings:
% .filtParams       filtParamHandler-object with settings for filter (more info: >>doc filtParamHandler)
% .reconParams      reconParamHandler-object with settings for reconstruction (more info: >>doc reconParamHandler)
% .visParams        visParamHandler-object with settings for visualization (more info: >>doc visParamHandler)
%
%Computation Times:
% .compTimeFilt     computation time for last execution of method "filter" [s]
% .compTimeRecon    computation time for last execution of method "reconstruct" [s] 
% .compTimeVis      computation time for last execution of method "visualize" [s] 
% .printCompTime    boolean deciding whether a computation time is printed on each call
%
%RX Data Size Properties (read only):
% .Na               number of transducer elements in .elemPos (Receiving transducer)
% .NaTX             number of transducer elements of transmitting probe (equals .Na unless .elemPosTx is defined)
% .Nacq             number of acquisitions (transmit events) in .TXPulseDelays
%
%Image Axes Related Properties (read only): 
% .Da               lateral width of the array [m] (if x-direction is lateral direction)
% .da               distance between elements in .elemPos [m]. If it varies by > 1 micro meter, calling .da causes an error
% .drs              spatial pulse-echo spacing according to temporal sampling frequency [m] (drs = c0/fs/2)
% .drc              spatial pulse-echo spacing according to temporal center frequency [m] (drs = c0/fc/2)
% 
%METHODS:
%
%Call:
%>> doc Recon2D.METHODNAME 
%for more information about the individual methods!
%
%Constructor and Copy:
% obj          = Recon2D(fs, elemPos, TXPulseDelays, c0, fc , profile)
% obj2         = obj.clone
%
%Data Processing: 
% [RXData, dt] = obj.filter(obj, RXData) 
% imgData      = obj.reconstruct(obj, RXData, dt, Z, X) 
% imgData      = obj.visualize(obj, imgData, Z, X)   
% imgData      = obj.run(obj, RXData, Z, X)
%
%Reconstructions:
% imgData      = obj.reconPWCuda(RXData_real, dt, z, x, recompileMex);
% imgData      = obj.reconPWCudaKSpace(RXData_real, dt, z, x, recompileMex);
% imgData      = obj.reconPBP(RXData, c, dz, dx);
% imgData      = obj.reconstructFixTXDelay(RXData, dt, Z, X, fixTXDelay)
% imgData      = obj.reconstructPlaneWaveTXRX(RXData, Z, X, dt, anglesTX, pivotPosTX, pivotDelayTX)
%
%Simulations:
% RXData       = obj.simulate(scatPos,nPwr,Nt,upsmplFct)
% RXData       = obj.simulateSerial(scatPos,nPwr,Nt)
% RXData       = obj.simulateFullTX(scatPos,nPwr,Nt,upsmplFct)
% RXData       = obj.simulateCUDA(scatPos,nPwr,Nt,simProps)
% scatPos      = obj.convertImage2scatPos(image,z,x,Nscat,params)
% RXData       = obj.simulateKWave(z,x,c,rho,simArgs)
% RXData       = obj.simulateKWave2(z,x,c,rho, plotSim,implMethod,sourceType)
%
%Advanced Data Processing:
% [r,x0,theta,err,maxElemErr] = obj.localizeArrayRadon(RXData, apertureID, settings)
%
%Parameter Processing:
%                obj.applyProfile(obj, profile)
%                obj.killFigs
%                obj.killDebugFigs
%                obj.setDebugFigs(hFilt,hRecon)
%
%Helpers:
%                obj.checkConsistency(RXData)
% waveField    = obj.plotWaveField(obj, z, x, t, i_acq)
%                obj.transformAperture(obj, T, apertureID)
%                obj.plotProbe()
% [z,x]        = obj.suggestAxes(spw_z, spw_x, Nz, Nx)
%                obj.plotTXDelays(RXData, i_acq, figHandle)
% TXDelays     = obj.setTXDelaysPlaneWave(angle, pivotDelay, pivotPos, normalAngle)
% TXDelays     = obj.setTXDelaysFocussedWave(focusPos)
%
%Acquisition Identification:
% arrayTypeStr = obj.getArrayType(apertureID, tol)
% isPW         = obj.isPlaneWave(tol)
% isPA         = obj.isPA()
% theta_tx     = obj.getPWAngles(tol)
%
%Converters:
% [Z,X]        = obj.getPolarGrids(obj,theta,rho,centerPos)
% [RXDataSynth, TXDelaysSynth] = obj.applySyntheticTXFocus(RXData, focPos)
% [RXDataSynth, TXDelaysSynth] = obj.applySyntheticTXPlaneWave(RXData, theta, i_pivot)
% virtualRXData = obj.virtualApertureRX(RXData, virtualRXElemPos, dt)
%
%AUTHOR: 
% Hans-Martin Schwab (hans-martin.schwab@web.de)

%TODO:
% - add a makeAxis method!!
% - memory save mode (only loops)
% - allow for cell array of filter types and params?
% - procDevice property sinnvoll? SCHNELLERER CPU CODE MOEGLICH ALS MIT ARRAYFUN? -> procDevice = 'CPU'; % string determining where data is processed (options:'CPU','GPU','CPUForced','GPUForced')
% - allow not to require all inputs in constructor and use dependent 'isReady' to verify ->  isReady %defines if all information required for the reconstruction is provided
% -complete commented methods
% -allow for stored geometry grids
% -complete getters, setters and privatize properties
% -make some methods static
% -check if it is faster to transfer all variables to GPU variables
% -write a method and a paramhandler-class-property for determinig T
% -setter that turns plotCompTime into 3 element array if scalar (each SP step asks forindividual element)
% -add element width to simulation by weighting signals according to
% elem-pixel angle 
% -if only one TX element exists, it will alwyas be interpreted as PA mode!
%
%TODO filter:
% - implement Nt_pad calculation for FD domain windows
% - improve antialias parametrization
% -write setters and getters of paramHandler
% -up/downmixing for IQDemod
% -allow for fullApertrueApodization
% -ADD FUNCTIONALITY FOR FIX COEFFS IN FILTPPARAMHANDLER?
% -add Bandlimit Parameter for visualization!
% -downsampling will be important to reduce GPU data transfer
% -reset funciton does not make sens ein this implementation
% -are the attenuation compensation values usually two-way??? Also: conversion correct (it's now .5@6dB)?
%
%TODO reconstruct:
% -write setters and getters of paramHandler
% -depthindependend apod with N aperture for faster comp?
% -adapt continous weight reconstruction
% -add memory save reconstruciton function?
% -also a fast version for CPU or a seperate one for stored geometry???
% -ADD FUNCTIONALITY FOR FIX DELAYS AND WEIGHTS in paramHandler?
% -wrapper for insoniWieght parameters (depending on pitch)
% -allow for directly setting default insoniWeigths in paramHandler
% -add option for empty cutoff angle to avoid computation of angle distance in recon
% -add property: .interpType {mustBeMember(interpType,{'nextNeighbor','linear'})} = 'nextNeighbor'; % interpolation type for RXData readout at delay times 
% -add property: .frameWeights = 1; % scalar weights for each frame (e.g. to scale contribution of TX angle)
% -elemDist not relaized! other parameterization!
% -write a reset funciton
%
%TODO visualize:
% -write setters and getters of paramHandler
% -is always pcolor just as fast?
% -write a reset funciton

classdef Recon2D < matlab.mixin.Copyable
%% PROPERTIES:   
    
properties
    
    % required recon properties:
    fs % sampling frequency [1/s]
    elemPos % Nx2 array with centerPositions of the elements [m] in order [zPos,xPos]
    TXPulseDelays % (NaTX x Nacq)-matrix with delays of the TX pulse peak for each element [s]
    c0 % scalar speed of sound [m/s]
    % optional recon properties:
    fc % center frequency [1/s] (not mandatory but required for some filters)
    normalAngleRX = 0; % .Na-element vector with normal angles of the elements (or scalar if equal)
    elemPosTX = [];% Nx2 array with centerPositions of the TX elements (if different) [m] in order [zPos,xPos] 
    % more optional properties:
    elemWidth % scalar containing the width of each transducer element [m]
    bandWidth % fractional bandwidth relative to .fc [1] (decimal, not in % !)
    
    % settings:
    filtParams % filtParamHandler-object with settings for filter (more info: >>doc filtParamHandler)
    reconParams % reconParamHandler-object with settings for reconstruction (more info: >>doc reconParamHandler)
    visParams % visParamHandler-object with settings for visualization  (more info: >>doc visParamHandler)

    % computation times:
    compTimeFilt  = nan; % computation time for last execution of method "filter" [s]
    compTimeRecon = nan; % computation time for last execution of method "reconstruct" [s] 
    compTimeVis   = nan; % computation time for last execution of method "visualize" [s] 
    printCompTime = false; % boolean deciding whether a computation time is printed on each call
    
    % other properties:
    gpuDevice % gpuDevice object if NVIDIA GPU can be found, empty otherwise
    freeParams = struct; % struct with freely chosable additional information
    
    %    %optionally storable geometry quantities:
    %METHODS:
    %.saveGeometry(Z, X, Nt)
    %.clearGeometry (must be called whenever a relevant param is changed?)
    %PROPERTIES (class geometryHandler?):
    %.geometryIsSaved % (DEPENDENT) boolean indicating if weights and delays are stored
    %.geometryDelays  % <- in combination with method ".saveGeometry(Z,X)"
    %.geometryApodWeights % <- in combination with method ".saveGeometry(Z,X)"
    %.geometryInsoniWeights  % <- in combination with method ".saveGeometry(Z,X)"

end

properties (Dependent, SetAccess = private)
    
    % data size properties:
    Na    % number of transducer elements in .elemPos (Receiving transducer)
    NaTX  % number of transducer elements of transmitting probe (equals .Na unless .elemPosTX is defined)
    Nacq  % number of acquisitions (transmit events) in .TXPulseDelays
    
    % axes related properties:
    Da   % lateral width of the array [m] (if x-direction is lateral direction)
    da   % distance between elements in .elemPos [m], or 'noEquidistantSpacing' if the difference in distance in > 1 micro meter
    drs  % spatial pulse-echo spacing according to temporal sampling frequency [m] (drs = c0/fs/2)
    drc  % spatial pulse-echo spacing according to temporal center frequency [m] (drs = c0/fc/2)
    
end

methods  
    
%% CONSTRUCTOR:

    function obj = Recon2D(fs, elemPos, TXPulseDelays, c0, fc, normalAngleRX, profile)
        %DASRECONSTRCUTOR delay-and-sum reconstruction class for US imaging 
        %for various configurations
        %   profile: 'fast','accurate','memorySave' 
        narginchk(0,7);        
        % check input:
        if nargin<1, fs = 4*8e6; 
         elseif isempty(fs), fs = 4*8e6;  end        
        if nargin<2, elemPos = [zeros(128,1),(-63.5:63.5)'*.3e-3]; 
         elseif isempty(elemPos), elemPos = [zeros(128,1),(-63.5:63.5)'*.3e-3]; end        
        if nargin<3, TXPulseDelays = zeros(size(elemPos,1),1); 
         elseif isempty(TXPulseDelays), TXPulseDelays = zeros(size(elemPos,1),1); end        
        if nargin<4, c0 = 1540;
         elseif isempty(c0), c0 = 1540;  end
        if nargin<5, fc = fs/4; 
         elseif isempty(fc), fc = fs/4; end
        if nargin<6, normalAngleRX = 0;
         elseif isempty(normalAngleRX), normalAngleRX = 0; end             
        if nargin<7, profile = 'fast'; end    
        % parse input:
        obj.fs = fs;
        obj.elemPos = elemPos;
        obj.TXPulseDelays = TXPulseDelays;        
        obj.c0 = c0;
        obj.fc = fc;
        obj.normalAngleRX = normalAngleRX;
        % initialize parameter handlers:
        obj.applyProfile(profile);  
        % connect GPU if parallel toolbox is found and NVIDIA gpu is found:
        if exist('gpuDeviceCount','file')==2 
            if gpuDeviceCount>0, obj.gpuDevice = gpuDevice; end 
        end
    end

%% METHODS - DATA PROCESSING:

    % filter data (definition in external file):
    [RXData, dt] = filter(obj, RXData);

    % reconstrcut data (definition in external file):    
    imgData = reconstruct(obj, RXData, Z, X, dt);

    % reconstrcut data (definition in external file):    
    imgData = reconstructVariousDelays(obj, RXData, Z, X, dt, taus); 

    % reconstrcut a coherence map (semblance) (definition in external file):    
    semblanceMap = reconstructSemblanceMap(obj, RXData, Z, X, dt, Npix, c_vals, applyAcqCompounding)   
    
    % reconstruct data using pre-computed arrival times (definition in
    % external file):
    imgData = reconstructPreCompRXDelay(obj, RXData, dt, Z, X, RXdelays, RXdelays_transmit);
    
    % reconstruct data using mex implementation:
    imgData = reconMex(obj,RXData_filt, z,x,dtu,recompileMex)
    
    % visualize data (definition in external file):    
    imgData = visualize(obj, imgData, Z, X); 
    
    % run all three steps at once:
    function imgData = run(obj, RXData, Z, X, bypassVis)
       %RUN Applies all three processing steps (filter, reconstruct, visualize)
       %consecutively.
       %INPUT:
       % - RXData: 2D or 3D array with receive data with temporal samples
       %           of all elements in columns and different transmission
       %           events in the third dimension, which will be compounded.
       % - Z: 1D array or 2D grid with z-positions of the image pixels [m]
       % - X: 1D array or 2D grid with x-positions of the image pixels [m]
       % - bypassVis: (optional) Boolean that allows to visualize log-
       %              compressed envelope data but to return RF data (or IQ 
       %              if chosen as filter output). Default value: false
       %OUTPUT:
       % - imgData: compounded image dataset on grid defined by z and x,
       %            equals output of .visualize (or of .reconstruct if 
       %            "bypassVis" is true.
       %AUTHOR: 
       %hans-martin.schwab@web.de (May2019)
       narginchk(4,5);
       % input check:
       if nargin<5, bypassVis = false; end % return vis output by default
       % filter, reconstruct, visualize
       [imgData, dtu] = obj.filter(RXData);
       imgData        = obj.reconstruct(imgData, Z, X, dtu);
       if bypassVis
            obj.visualize(imgData, Z, X); 
       else
            imgData   = obj.visualize(imgData, Z, X);            
       end
       % print computation time if desired:
       if obj.printCompTime
           fprintf('-> total computation time:      < %.3fs > \n',...
                   obj.compTimeFilt+obj.compTimeRecon+obj.compTimeVis);
       end
    end

%% METHODS - ADVANCED DATA PROCESSING:

    % plane wave CUDA DAS Recon:
    imgData = reconPWCuda(obj, RXData_compl, dt, z, x, params, recompileMex);
    
    % plane wave CUDA k-space Recon:
    [imgData, z ,x] = reconPWCudaKSpace(obj, RXData_real, dt, params, recompileMex);

    % ReconWithFixTXDelay:
    imgData = reconstructFixTXDelay(obj, RXData, dt, Z, X, fixTXDelay);

    % plane wave DAS recon with plane RX angles:
    [imgData, RXAnglData] = reconstructPlaneWaveTXRX(obj, RXData, Z, X, dt, anglesTX, pivotDelayTX, pivotPosTX, anglesRX)

    % reconstruction with aberration correction using PBP:
    [imgData, fieldRX, fieldTX] = reconPBP(obj, RXData, c, dz, dx);
    
    % plane wave linear array localization for dual probe imaging:
    [r,x0,theta,err,maxElemErr] = localizeArrayRadon(obj, RXData, apertureID, settings)
    
%% METHODS - DATA MANAGMENT: 

    %     function deleteGeometry(obj, fieldsStringCellArray)
    %         % delete optionally storable geometry quantities.
    %         % 
    %     end

    %     function saveGeometry(obj, Z, X)
    %         % store storable quantities
    %         % 
    %     end
    
    % check if sizes of properties are consitent:
    function checkConsistency(obj, RXData)
       %CHECKCONSISTENCY checks if the sizes of the properties .elemPos, 
       % .TXDelays and the (optional) input argument "RXData" are consitent
       % meaning refer to the same number of elements, and the same number
       % of transmit events. The function does not return anything but
       % raises an error in case of inconsitency.
       %AUTHOR: hans-martin.schwab@web.de (June20)
       narginchk(2,3);
       % check consitency of TX elements: 
       if ~isempty(obj.elemPosTX)
            % only check number of TX delays if no PA mode:
            if size(obj.TXPulseDelays,1)~=1
            assert(size(obj.elemPosTX,1)==size(obj.TXPulseDelays,1), ...
                   ['Inconsitency detected: The number of transmit elements in .elemPosTX is ',...
                    num2str(size(obj.elemPosTX,1)),', but the number of transmit elements in .TXPulseDelays is ',...
                    num2str(size(obj.TXPulseDelays,1)),'!']); 
            end
       else
            % only check number of TX delays if no PA mode:
            if size(obj.TXPulseDelays,1)~=1           
            assert(size(obj.elemPos,1)==size(obj.TXPulseDelays,1), ...
                   ['Inconsitency detected: The number of transmit elements in .elemPos is ',...
                    num2str(size(obj.elemPos,1)),', but the number of transmit elements in .TXPulseDelays is ',...
                    num2str(size(obj.TXPulseDelays,1)),' (and no extra .elemPosTX is defined)!']);         
            end
       end
       if nargin==2
            % check consitency of RX elements: 
            assert(size(obj.elemPos,1)==size(RXData,2), ...
                   ['Inconsitency detected: The number of receive elements in .elemPos is ',...
                   num2str(size(obj.elemPos,1)),', but the input dataset RXData has ',...
                   num2str(size(RXData,2)),' columns!']);
            % check consitency of acquisitions:   
            assert(size(obj.TXPulseDelays,2)==size(RXData,3), ...
                   ['Inconsitency detected: The number of acquisitions in .TXPulseDelays is ',...
                   num2str(size(obj.TXPulseDelays,2)),', but the input dataset RXData has ',...
                   num2str(size(RXData,3)),' elements in the 3rd dimension!']);               
       end
    end

%% METHODS - SETTINGS MANAGEMENT:

    function applyProfile(obj, profile)
        % sets the processing settings to a predifined profile
        % options: 'fast','accurate','memorySave' 
        obj.filtParams  = filtParamHandler(profile);
        obj.reconParams = reconParamHandler(profile);
        obj.visParams   = visParamHandler(profile);  
    end
    
    % clone:
    function copiedObj = clone(obj)
        %CLONE clone the object to get a new object with same properties.
        narginchk(1,1);
        % copy class object:
        copiedObj = obj.copy;
        % copy all objects used in this class individually:
        copiedObj.filtParams = obj.filtParams.copy;
        copiedObj.reconParams = obj.reconParams.copy;
        copiedObj.visParams = obj.visParams.copy;    
    end

    % deactivate all figure handles:
    function killFigs(obj)
        %KILLFIGS deactivate all figure handles for filter, recon and visualization
       obj.killDebugFigs;
       obj.visParams.figHandle = [];
    end
    
    % deactivate figures:
    function killDebugFigs(obj)
       %KILLDEBUGFIGS deactivate figures for filter and recon 
       obj.filtParams.figHandle = [];
       obj.reconParams.figHandle = [];
    end

    % activate figures for filter and recon with handles or figure indeicees:
    function setDebugFigs(obj,hFilt,hRecon)
       %SETDEBUGFIGS activate figures for filter and recon with handles or figure indeicees
       % INPUT:
       % - hFilt: figure handle or figure ID for filter visualization
       % - hRecon: figure handle or figure ID for recon visualization       
       if nargin<2, hFilt = figure; end
       if nargin<3, hRecon = figure; end
       obj.filtParams.figHandle = figure(hFilt);
       obj.reconParams.figHandle = figure(hRecon);
    end
     
    % check if gpu is connected:
    function val = checkForGPU(obj)
        %CHECKFORGPU checks if an NVIDIA GPU has been found and can be
        %found in the property .gpuDevice.
        % If the method is called with and output argument the output is 
        % boolean.
        % If it is called without an output argument it does nothing if
        % true but raises an error if wrong.
        %AUTHOR: hans-martin.schwab@web.de (Jan2022)
        val = true;
        if isempty(obj.gpuDevice), val = false; end
        if nargout==0, assert(val, ['This method cannot be run, because no ',...
        	'suitable GPU was found, or the Parallel Computing Toolbox is not installed! ']); 
        end
    end
    
%% METHODS - HELPERS: 

    %     function dr = getNyqusitSpacing(obj, fc, bw) % ACHTUNG: BEI RF BESSER fs ANGEBEN? -> HERAUSFINDEN!
    % 
    %     end

    %     %static:
    %     function sensiCurve = getElemDirectivityCurve(obj, elemWidth, fc, theta)
    % 
    %     end

    % plot transducers into image:
    function plotProbe(obj, color, plotLegend)
        % PLOTPROBE plots the position the transducer elements into the current
        % figure (on top of previously visualized content in that figure.
        % If transmit and receive aperture differ, both are plottet in
        % different colors. If the normal angles of the receive aperture are
        % defined (in the property .normalAngleRX), the direction is
        % also visualized as one line per receive element. The optional input
        % argument "color" is a color in MATLAB color specification format 
        % (RGB triplet, shirt name string or long name string) or a two 
        % element cell array {TX,RX} for different TX and RX probe positions.
        % Another optional input argument "plotLegend" is a boolean
        % defining whether a legend should be shown, which is false by
        % default.
        %AUTHOR: hans-martin.schwab@web.de
        narginchk(1,3);
        % input check:
        if nargin < 2, color={'m','g'}; elseif isempty(color), color={'m','g'}; 
         elseif ~iscell(color), color = {color,'g'}; end
        if nargin < 3, plotLegend = false; end
        % plot transducer:
        hold all; scatter(obj.elemPos(:,2)*1e3,obj.elemPos(:,1)*1e3, [color{2},'.']);  
        legendText{1} = 'transducer elements';
        % plot TX transducer (if different):
        if ~isempty(obj.elemPosTX)
            scatter(obj.elemPosTX(:,2)*1e3,obj.elemPosTX(:,1)*1e3,[color{1},'.']);
            legendText{end}   = 'transducer elements (RX)';
            legendText{end+1} = 'transducer elements (TX)';
        end
        if isscalar(obj.normalAngleRX)
           theta = obj.normalAngleRX*ones(obj.Na,1);
        else
           theta = obj.normalAngleRX; 
        end
        % set length of the line plot as 5x element distance:
        lineLength = 4*sqrt(sum(diff(obj.elemPos(1:2,:),1,1).^2))*1e3;
        % plot normal angles of elements:
        plot(obj.elemPos(:,2)'*1e3 + lineLength*[zeros(obj.Na,1),sind(theta)]',...
             obj.elemPos(:,1)'*1e3 + lineLength*[zeros(obj.Na,1),cosd(theta)]',...
             color{2},'HandleVisibility','off');       
        if plotLegend, legend(legendText,'Location','south'); end
        hold off;
    end
    
    function transformAperture(obj, T, apertureID)
        %TRANSFORMAPERTURE changes the probe's location and orientation by
        % applying a geometric tranfrom to all elements and by rotating the
        % element's normal angles.
        %INPUT:
        % - T: either: 3x3 geometric transform matrix of the shape: 
        %              [ cosd(theta), sind(theta), 0;
        %               -sind(theta), cosd(theta), 0;
        %                       tz_p,        tx_p, 1] 
        %              such that it can be multiplied with a position 
        %              vector of the shape: [z,x,1]*T 
        %      or: cell array with 3 or 5 elements parameterizing T as:
        %          {tz,tx,theta, tz_p,tx_p}, where tz and tx determine the 
        %          translation before rotation and the optional parameters 
        %          tz_p and tx_p (p=prime) determine a translation after rotation.
        % - aptertureID: string that defines, which aperture should be
        %                transformed , e.g. for multiprobe imaging. Options: 
        %                -'both'(default) tranforms RX and TX aperture if defined 
        %                -'RX' transforms only RX aperture and keeps TX aperture
        %                      (if  no TX aperture was defined, it is
        %                      created as as current RX aperture before transform)        
        %                -'TX' transforms only TX aperture and keeps RX aperture
        %                      (if  no TX aperture was defined, it is
        %                      created as transform of current RX aperture)
        %AUTHOR: Hans-Martin Schwab
        narginchk(2,3);
        % input checking:
        if nargin<3
           apertureID = 'both'; 
        end
        % create transformation matrix if given in parameters:
        if iscell(T)
            assert(length(T)==3 | length(T)==5,...
                   'If input argument T is a cell array it must have 3 or 5 elements: {tz,tx,theta, tz_p,tx_p}')
            % set default translations after rotation to zero if not defined:
            if length(T)==3, T{4}=0; T{5}=0; end
            % create transform matrix:
            T = [ cosd(T{3}), sind(T{3}), 0;...
                 -sind(T{3}), cosd(T{3}), 0;...
                 T{1}*cosd(T{3})-T{2}*sind(T{3})+T{4},...
                 T{2}*cosd(T{3})+T{1}*sind(T{3})+T{5}, 1];
        end
        % check T:
        assert(all(size(T)==[3,3]), 'Input argument T must be a 3x3 matrix!');
        % check aperture ID:
        assert(any(strcmp(apertureID,{'RX','TX','both'})),...
               "Input argument 'apertureID' must be 'RX','TX' or 'both'!");
        % transform RX aperture if desired:
        if any(strcmp(apertureID, {'RX','both'}))
                %create a TX aperture if undefined before transforming only RX aperture:
                if strcmp(apertureID,'RX') && isempty(obj.elemPosTX)
                    obj.elemPosTX = obj.elemPos; end
                %check if RX aperture was defined:
                assert(~isempty(obj.elemPos),'Property .elemPos is empty and cannot be transformed');
                % transform RX aperture:
                newElemPos = [obj.elemPos, ones(obj.Na,1)] * T;
                % set RX aperture to transformed aperture:
                obj.elemPos = newElemPos(:,1:2);
                % rotate normal angle RX by transform angle:
                rotNormalVect = ...
                        [cosd(obj.normalAngleRX(:)),...
                         sind(obj.normalAngleRX(:))] * T(1:2,1:2);
                obj.normalAngleRX = atan2d(rotNormalVect(:,2),rotNormalVect(:,1));
        end
        % transform TX aperture if desired:
        if any(strcmp(apertureID, {'TX','both'}))
                % if not TX aperture is defined, transform RX apert. to TX apert.:
                if isempty(obj.elemPosTX)
                   % transform RX aperture:
                   newElemPos = [obj.elemPos, ones(obj.Na,1)] * T;
                else
                    % transform TX aperture:
                    newElemPos = [obj.elemPosTX, ones(obj.Na,1)] * T;                    
                end
                % set TX aperture to transformed aperture:
                obj.elemPosTX = newElemPos(:,1:2);
        end
    end


    % apply a serial time-of-flight simulation (slow):
    RXData = simulateSerial(obj, scatPos, nPwr, Nt);
    
    % a quicker, slightly less accurate version of the simulateSerial method:
    RXData = simulate(obj, scatPos, nPwr, Nt, upsmplFct);
    
    % a slower, more accurate version of the simulateSerial that takes the 
    % full TX wavefront into account:
    RXData = simulateFullTX(obj, scatPos, nPwr, Nt, upsmplFct);
       
    % a fast, accurate and more model based CUDA TOF simulation:
    RXData = simulateCUDA(obj, scatPos, nPwr, Nt, simProps);
    
    % a pixel based pseudospectral simulation (k-wave):
    RXData = simulateKWave(obj,z,x,c,rho,simArgs)
    
    % a pixel based pseudospectral simulation (k-wave) using the PULS/e k-wave toolbox:
    RXData = simulateKWave2(obj,z,x,c,rho, plotSim,implMethod,sourceType)
    
    % static:
    % read image data and create scatter Position data
    scatPos = convertImage2scatPos(obj,image,z,x,Nscat,params)
    
    % plot the wavefield accoring to geometry and TXDelays:
    function waveField = plotWaveField(obj, z, x, t, i_acq)
        %PLOTWAVEFIELD plots the wavefield that is emitted according to the 
        % element positions and the transmission delays on a regular grid
        % defined by the axes z and x on time steps t and returns the
        % wavefield in a 3D array
        %INPUTS:
        % - z: 1D array with z-axis on which to display wavefield [m]
        % - x: 1D array with z-axis on which to display wavefield [m] 
        % - t: 1D array with time steps at which to display wavefield [s]
        % - i_acq: interger with index of TX event which shall be plotted
        %          (column index of .TXPulseDelays)
        %OUTPUT:
        % - waveField: 3D array with wavefields on z,x axes for each time
        %              step in t of size [Nz,Nx,Nt]
        %TODO: allow z and x to be grids
        narginchk(4,5);
        % input checking:
        assert(isvector(x)&&isvector(z)&&isvector(t), 'Input arguments z, x and t must be vectors!');
        if nargin<5, i_acq=1; end
        assert(i_acq>=0&&i_acq<=obj.Nacq, 'input argument i_acq is not in the range [1,.Nacq]!');
        % precomputation:
        if isempty(obj.elemPosTX), TXElemPos = obj.elemPos; else, TXElemPos = obj.elemPosTX; end
        zElem = permute(TXElemPos(:,1),[3,2,1]);
        xElem = permute(TXElemPos(:,2),[3,2,1]);
        tauPulse = permute(obj.TXPulseDelays(:,i_acq),[3,2,1]);
        % wave field computation for each time step:
        if nargout>0, waveField = zeros(length(z),length(x),length(t)); end
        for i_t = 1:length(t)
            % compute input time minus time from element to pixel:
            pixelDelays = t(i_t) - tauPulse - sqrt((z(:)-zElem).^2+(x(:)'-xElem).^2)/obj.c0;
            % generate wavefield by superposing waves from all elements:
            waveField_i = sum(gauspuls(pixelDelays, obj.fc),3);
            % display wavfield:
            imagesc(x*1e3,z*1e3, waveField_i); axis image; colormap cool; drawnow;
%             obj.plotProbe; xlabel('x [mm]'); ylabel('z [mm]'); 
            % save wavefield to output if desired:
            if nargout>0, waveField(:,:,i_t) = waveField_i; end
        end
    end
    
    %suggest axes according to center fequency and probe dimesions:
    function [z,x] = suggestAxes(obj, dz, dx, Nz, Nx)    
        %SUGGESTAXES returns axis arrays based on the (optional inputs)
        %spacings dz and dx and the optional lengths Nz and Nx. Axes are 
        %centered with the probe center in lateral direction and start at
        %the minimum element depth in axial direction.
        %By default, dz and dx are 2 and 1 samples per center wavelength
        %and Nz and Nx are chosen, such that both axes match the probe width.
        %(Note that .elemPos, .fc and .c0 must be defined to run this method!) 
        %AUTHOR: hans-martin.schwab@web.de (Dec19)
        %TODO: include arbitrary cases where aperture is rotated!
        narginchk(1,5);
        % get properties depending on probe dimensions:
        xc = min(obj.elemPos(:,2)) + obj.Da/2; % lateral probe center position
        z0 = min(obj.elemPos(:,1)); % minimum axial element position
        % set input defaults:
        if nargin<2, dz = 1*obj.drc; elseif isempty(dz), dz = 1*obj.drc; end % 2 samples per center wavelength   
        if nargin<3, dx = 2*obj.drc; elseif isempty(dx), dx = 2*obj.drc;  end % 1 sample per center wavelength    
        if nargin<4, Nz = round(obj.Da/dz)+1; end % grid width equals probe width
        if nargin<5, Nx = round(obj.Da/dx)+1; end % grid depth equals probe width
        % create axes:
        z = z0 + (0:Nz-1)'*dz;
        x = xc + (-floor(Nx/2):floor((Nx-1)/2))*dx;
    end
    
    % plot TXDelays into RXData for verification:
    function plotTXDelays(obj, RXData, i_acq, figHandle)
        %PLOTTXDELAYS plots the TXDelays into the RXData set for verification
        %INPUTS:
        % - RXData:     (optional) 3D array with Receive data with dimension [Nt,.Na,.Nacq] 
        % - i_acq:      (optional) scalar integer with index of acquisition 
        %                in RXData that shall be displayed (default: 1)
        % - figHandle: (optional) handle to the figure that shows the plot,
        %               creates new figure by default
        % AUTHOR: hans-martin.schwab@web.de
        narginchk(1,4);
        % input parsing:
        if nargin<2, RXData = []; end
        if nargin<3, i_acq = 1; end % use first acquistion by default
        if nargin<4, figHandle = figure; end % create new figure if no handle is specified
        assert(isscalar(i_acq),'Input argument "i_acq" must be scalar!');
        assert(i_acq>=1 & i_acq<=obj.Nacq & mod(i_acq,1)==0,...
               'Input argument "i_acq" must be an integer between 1 and .Nacq!');
        % create figure:
        figure(figHandle);
        % plot RX Data if provided:
        if ~isempty(RXData)
            imagesc(1:obj.Na, (0:size(RXData,1)-1)/obj.fs*1e6, single(real(RXData(:,:,i_acq)))); colormap gray;
            ylim([min([0;obj.TXPulseDelays(:,i_acq)]), (size(RXData,1)-1)/obj.fs]*1e6);
        end
        % plot TX delays:
        hold all; scatter(1:obj.Na, obj.TXPulseDelays(:,i_acq)*1e6,'g.'); hold off;
        xlabel('channel'); ylabel('time [\mus]');
    end

%% METHODS - ACQUISITION IDENTIFICATION:

    function arrayTypeStr = getArrayType(obj, apertureID, tol)
        %ISLINARRAY returns type of array and returns it as a string
        % (options: "linear", "curved", "other"). The optional input
        % argument apertureID defines which aperture is checked, if transmit
        % and receive apertures differ (property .elemPosTX is not empty) 
        % (options: "TX","RX"[default]). The optional input argument
        % tolerance "tol" [m] defines the tolerance by which an element
        % distance might differ from the exact position for that type,
        % the default value is 1e-6 (1 micro meter).
        %AUTHOR: hans-martin.schwab@web.de (Jan2020)
        narginchk(1,3);
        % set defaults:
        if nargin < 2, apertureID = "RX"; end
        if nargin < 3, tol = 1e-6; end
        % write element positions of desired apterture into "pos":
        assert(any(strcmp(apertureID,{'RX','TX'})),'Input argument apertureID must be "RX" or "TX"');        
        if strcmp(apertureID,'RX')||isempty(obj.elemPosTX), pos = obj.elemPos; else; pos = obj.elemPosTX; end
        % identify as linear array, if the distances between neighboring 
        % elements in both in x- and z- direction are all the same:
        if all(abs(diff(pos(:,1),2,1)) < tol) && all(abs(diff(pos(:,2),2,1)) < tol)
            arrayTypeStr = 'linear';
        % identify as curved array, if the Eucledean distances between neighboring 
        % elements are all the same and the angle increase between two neighboring 
        % connection lines of two elements are all the same (there, the tolerance 
        % is defined as the angle difference, when the next element is off by "tol"):            
        elseif  all(abs(diff(sqrt(diff(pos(:,1)).^2 + diff(pos(:,2)).^2))) < tol)...
              && all(diff(atan2(diff(pos(:,1)),diff(pos(:,1)))) ...
                     < atan2(tol,sqrt(diff(pos(1:2,1)).^2 + diff(pos(1:2,2)))))             
            arrayTypeStr = 'curved';
        else
            arrayTypeStr = 'other';
        end
    end
    
    % check if acquisitoin is photoacoustic (only one TX delay per acquisition)
    function isPA = isPA(obj)
       if size(obj.TXPulseDelays,1)==1, isPA = true; else, isPA = false; end 
    end
    
    function isPW = isPlaneWave(obj, tol)
        %ISPLANEWAVE returns a boolean that is true if the transmit wave
        % defined by the property .TXPulseDelays is a plane wave. The
        % optional input argument "tol" defines the tolarance between the
        % change of transmit delay differences over the array to still
        % count as plane wave in seconds (default: 10e-9). For multiple 
        % transmissions (.Nacq>1), all transmissions must be plane waves. 
        % Note that this method only works for linear arrays!
        %AUTHOR: hans-martin.schwab@web.de (Jan2020)
        narginchk(1,2);
        if nargin < 2, tol = 10e-9; end 
        assert(strcmp(obj.getArrayType('TX'),'linear'), ['The method .getPWAngles ',...
               'is only implemented for linear arrays but not linear array ',...
               'is defined in the property .elemPos / .elemPosTX!']);        
        % identify as plane wave if differences between transmission delays
        % increase linearly with element
        isPW = all(all(abs(diff(obj.TXPulseDelays,2,1)) < tol ));
    end
    
    function theta_TX = getPWAngles(obj)  
        %GETPWANGLES returns the plane wave TX angles based on the transmission 
        % delays in the property .TXPulseDelays as an array of angles per
        % acquisisiton (in case they describe a plane wave transmission). 
        % Note that this method only works for linear arrays!
        %AUHTOR: hans-martin.schwab@web.de (Jan2020)
        assert(strcmp(obj.getArrayType('TX'),'linear'), ['The method .getPWAngles ',...
               'is only implemented for linear arrays but not linear array ',...
               'is defined in the property .elemPos / .elemPosTX!']);
        assert(obj.isPlaneWave, ['No transmit angle could ',...
               'be found, because the delays in the property .TXPulseDelays ',...
               'do not refer to a plane wave!']);
        theta_TX = 90 - acosd((obj.TXPulseDelays(2,:) - obj.TXPulseDelays(1,:))...
                        *obj.c0 / obj.da);
    end
    
%% METHODS - CONVERTERS: 

    % static
    function [Z,X] = getPolarGrids(obj,theta,rho,centerPos)
        %GETPOLARGRIDS converts vectors of angles "theta" [deg] and radii 
        % "rho" [m] and an optional center position "centerPos" [z,x] 
        % (default: [0,0]) into carthesian grid matrixes Z and X [m] according to
        % the class geometry definitions (see drawing in >>doc Recon2D)
        %AUTHOR: hans-martin.schwab@web.de
        narginchk(3,4);
        % input check:
        if nargin < 4
           centerPos = [0,0]; 
        end    
        % get full grids:
        [Rho,Theta] = ndgrid(rho,theta/180*pi);
        % get carthesian coordinates of these gridpoints:
        [Z,X] = pol2cart(Theta,Rho);
        % translate to desired center:
        Z = Z + centerPos(1);
        X = X + centerPos(2);    
    end

    function TXDelays = setTXDelaysPlaneWave(obj, angle, pivotDelay, pivotPos, normalAngle)
        % SETTXDELAYSPLANEWAVE computes the transmit delays for the
        % current transmit aperture that describe a plane wave according 
        % to the parameters in the input arguments and either returns
        % them (if an output argument is defined) or sets the property
        % .TXPulseDelays to these delays (if no output arguments is defined)
        %
        %INPUT:
        % - angle: vector with transmit angles (direction wavefront is 
        %          travelling relative to input argument "normalAngle") 
        %          in [deg] (counter-clockwise)
        % - pivotDelay: (optional) time between acquisition start and 
        %               transmission of wave at pivot position defined in
        %               input argument "pivotPos" in [s] (default: 0)
        % - pivotPos: (optional) position of the pivot of the plane wave,
        %             either a scalar (pivot element index) or a two-element 
        %             vector (pivot position [z,x]), can also be a 
        %             (Nacq x 1/2)-vector/matrix, if pivot point changes
        %             with the angle (default: 1);
        % - normalAngle: (optional) normal angle of aperture, angle in
        %                global coordinate, relative to which the input
        %                argument "angle" is defined in [deg] (default: 0)
        %OUTPUT: 
        % - TXDelays: (optional) (Na x Nacq) matrix with transmit delays 
        %             for each acquisition and each angle (if output is not
        %             returned, the object's TXDelays will be overwritten)
        narginchk(2,5);
        % input check:
        if nargin<3
           pivotDelay = 0; 
        end
        if nargin<4
            pivotPos = 1;
        end
        if nargin<5
            normalAngle = 0;
        end
        % get positions of transmit aperture:
        if isempty(obj.elemPosTX)
           elementPos = obj.elemPos;
        else
           elementPos = obj.elemPosTX; 
        end        
        % convert pivot to position if it is an element index:
        if isscalar(pivotPos)
           assert(pivotPos>=1 & pivotPos<=obj.Na & mod(pivotPos,1)==0,...
               'If input argument "pivotPos" is scalar it must be an integer >=1 and <=.Na!') 
           pivotPos = elementPos(pivotPos,:);
        else
           assert(numel(pivotPos)==2,'Input argument "pivotPos" must be a scalar or a two-element vector!') 
        end
        % compute delays:
        absoluteAngle = normalAngle + angle(:)';
        TXDelays = pivotDelay(:)' + (elementPos(:,1)-pivotPos(1,1)')...
                                   *cosd(absoluteAngle)/obj.c0...
                                  + (elementPos(:,2)-pivotPos(1,2)')...
                                   *sind(absoluteAngle)/obj.c0;
        % if no output argument is defined, the object's delays are set:
        if nargout==0, obj.TXPulseDelays = TXDelays; end
    end

    % set transmit delays for a focussed wave emission:
    function TXDelays = setTXDelaysFocussedWave(obj, focusPos, subApertureTX)
        % SETTXDELAYFOCUSSEDWAVE computes the transmit delays for the
        % current transmit aperture that describe a focussed wave according 
        % to the parameters in the input arguments and either returns
        % them (if an output argument is defined) or sets the property
        % .TXPulseDelays to these delays (if no output arguments is defined)
        % .TXPulseDelays are focussed waves if the specified focus position 
        % is behind the aperture, the focus will be mirrored at the transducer.
        %
        %INPUT:
        % - focusPos: Npos-by-2 matrix with focus positions in [z,x] order
        % - subApertureTX: (optional) either:
        %                    - binary Na-by-Nfoci matrix defining the
        %                      active transmit elements for each focus position
        %                   or:
        %                    - scalar representing the F-Number (ratio between 
        %                    focus depth and aperture size. The subaperture
        %                    is calculated from this (Note that the property 
        %                    .normalAngleRX must be set correctly for this
        %                    option to work!).
        %                  (default: ones(Na,Nfoci))
        %
        %OUTPUT: 
        % - TXDelays: (optional) (Na x Nacq) matrix with transmit delays 
        %             for each acquisition and each angle (if output is not
        %             returned, the object's TXDelays will be overwritten)
        %AUTHOR: hans-martin.schwab@web.de
        narginchk(2,3);
        % input check:
        if nargin<3, subApertureTX = ones(obj.Na,size(focusPos,1)); end
        % create binary subAperture matrix if defined by scalar F-Number:
        if isscalar(subApertureTX)
            subApertureTX = abs(atan2(focusPos(:,2)'-obj.elemPos(:,2),...
                                      focusPos(:,1)'-obj.elemPos(:,1))...
                                - obj.normalAngleRX(:)/180*pi)...
                            <= atan(2*subApertureTX) ;
        end
        assert(size(focusPos,2)==2,...
               'The input argument "focusPos" must be a two-column array!')
        assert(all(size(subApertureTX)==[obj.Na,size(focusPos,1)]),...
               'The input argument subAperture must be binary a Na-by-Nfoci array!');
        % get positions of transmit aperture:
        if isempty(obj.elemPosTX), elemPosition = obj.elemPos; 
        else, elemPosition = obj.elemPosTX; end   
        % calculate delays:
        TXDelays = - sqrt( (elemPosition(:,1)-focusPos(:,1)').^2 ...
                          +(elemPosition(:,2)-focusPos(:,2)').^2) / obj.c0; 
        % apply subaperture masking:
        TXDelays(~subApertureTX) = inf;
        % correct for minimum delay:
        TXDelays = TXDelays - min(TXDelays,[],1);        
        % if no output argument is defined, the object's delays are set:
        if nargout==0, obj.TXPulseDelays = TXDelays; end
    end    
    
    % shift RX data to apply synthetic plane wave:
    [RXDataSynth, TXDelaysSynth] = applySyntheticTXPlaneWave(obj, RXData, theta, pivotPos)
    
    % shift RX data to apply synthetic TX foci:
    [RXDataSynth, TXDelaysSynth] = applySyntheticTXFocus(obj, RXData, focPos);
    
    % generate virtual receive data from arbitrary positions:
    virtualRXData = virtualApertureRX(obj, RXData, virtualRXElemPos, dt)
    
    % crate a circular array (.elemPos) for IVUS / curved arrays
    function [elemPos, normalAngleRX] = createCircularArray(obj, Na, radius, arcAngle, centerPos, isInward)
       %CREATECIRCULARARRAY creates element positions and normal angles
       %that describe (an arc of) a circular aperture. If no outputs are
       %defined, the element positions and normal angles are changed,
       %otherwise they remain uneffected!
       %USAGE:
       % [elemPos, normalAngleRX] = createCircularArray(obj, Na, radius, arcAngle, centerPos, isInward)
       %INPUT:
       % - Na:          scalar with number of transducer elements [1]
       % - radius:      scalar with radius of circular aperture [m]
       % - arcAngle:    scalar or 2 element array with full angle symmetric
       %                around 0 or start and stop angle of aperture [deg] (default: 360)
       % - centerPos:   2 element vector with center position [z,x] [m] (default: [0,0]) 
       % - isInward:    boolean detemrining if normal angle points inward (tomography)
       %                or outbound (curved array) (defaulkt: false)
       %OUTPUT:
       % - elemPos:       element position defined as in property .elemPos 
       %                  (property unchanged if output is defined)
       % - normalAngleRX: normal angles defined as in property .normalAngleRX
       %AUTHOR: hans-martin.schwab@web.de (2021)
        narginchk(3,6);
        % parse input:
        if nargin<4, arcAngle = 360; elseif isempty(arcAngle), arcAngle = 360; end 
        if nargin<5, centerPos = [0,0]; elseif isempty(centerPos), centerPos = [0,0]; end 
        if nargin<6, isInward = false; end
        if isscalar(arcAngle), arcAngle = arcAngle/2*[-1,1]; end
        % define center-to-element angles (avoid redundant positoin for 360 deg)!
        if abs(diff(arcAngle)) == 360
            normalAngleRX = arcAngle(1) + (0:Na-1)'/Na*diff(arcAngle);
        else
            normalAngleRX = linspace(arcAngle(1),arcAngle(2),Na)';
        end
        % get element postions
        elemPos = radius*[cosd(normalAngleRX),sind(normalAngleRX)] + centerPos;
        % rotate normal angles inward if desired:
        if isInward, normalAngleRX = mod(normalAngleRX+180, 360); end
        % only write to obj if no output is defined:
        if nargout<1, obj.elemPos = elemPos; obj.normalAngleRX = normalAngleRX; end
    end

%% SETTERS AND GETTERS:

    % Required Properties:
    %elemPos:
    function set.elemPos(obj,val)
        assert(isempty(val)||size(val,2)==2,...
               'the property ".elemPos" must be a two-column matrix [z,x]');
        obj.elemPos = val;
    end    
    
    % Optional Properties:
    %fc:
    function val = get.fc(obj)
        assert(~isempty(obj.fc),'the property .fc must be defined to proceed this operation!');
        val = obj.fc;
    end    
    %elemPosTX:
    function set.elemPosTX(obj,val)
        assert(isempty(val)||size(val,2)==2,...
               'the property ".elemPos" must be a two-column matrix [z,x]');
        obj.elemPosTX = val;
    end
    %TXPulseDelays:
    function set.TXPulseDelays(obj,val)
        assert(isempty(val)||size(val,1)==obj.Na||size(val,1)==1,...
               ['The property ".TXPulseDelays" have .Na rows (for US mode) ',...
                'or a single row (for PA mode)']);
        obj.TXPulseDelays = val;
    end
    
    % Dependent Properties:
    %Na:
    function val = get.Na(obj),   val=size(obj.elemPos,1); end 
    %Nacq:
    function val = get.Nacq(obj), val=size(obj.TXPulseDelays,2); end 
    %NaTX:
    function val = get.NaTX(obj)
        if isempty(obj.elemPosTX)
            val=size(obj.elemPos,1); 
        else
            val=size(obj.elemPosTX,1);            
        end
    end
    %drs:
    function val = get.drs(obj), val=obj.c0/obj.fs/2; end 
    %drc:
    function val = get.drc(obj), val=obj.c0/obj.fc/2; end 
    %Da: % lateral probe width 
    function val = get.Da(obj), val = max(obj.elemPos(:,2))-min(obj.elemPos(:,2)); end 
    %da:
    function val = get.da(obj)
        % check if at least two elements are defined!
        assert(obj.Na>1, 'There must be at least 2 elements in .elemPos to determine .da!');
        da_all = sqrt(diff(obj.elemPos(:,1)).^2 + diff(obj.elemPos(:,2)).^2);
        % check if distance between elements is the same 
        %(tolerance of 1 micro meter to avoid rounding errors etc.)
        assert(all(abs(diff(da_all)) < 1e-6), ['.da could not be determined,',...
               ' because the distance between elements varies within .elemPos by >1mu!']);
        % assign values of first distcance:
        val = da_all(1); 
    end    
    
end
    
end
