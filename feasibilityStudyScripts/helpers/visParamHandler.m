%RECONPARAMHANDLER Class to manage visualization parameters of DASReconstructor class
%
%PROPERTIES:
%.dynamicRange   = [];    % scalar positive value for in dB for logarithmic compression (if []: linear scale)
%.applyEnvelope  = true;  % boolean determining if envelope detection should be applied   
%.padFactHilbert = 1;     % factor by which signal is padded before Hilbert trafo
%.isEquidistGrid = true;  % boolean determining if the image grid is equidistant
%.figHandle      = 'last';% handle to figure showing the image (other options: [](no fig),'new','last')
%.figTitle       = '';    % title of figure (if .figHandle=[]: no effect)
%.supressDrawnow = false; % boolean deciding of figure is immediately plotted, set true to add more curves to figure;
%
%AUTHOR: 
% Hans-Martin Schwab (hans-martin.schwab@web.de)

classdef visParamHandler < matlab.mixin.Copyable
    %RECONPARAMHANDLER Class to manage visualization parameters of DASReconstructor class
    
properties
    
    dynamicRange = 60; % scalar positive value for in dB for logarithmic compression (if []: linear scale)
    applyEnvelope = true; % boolean determining if envelope detection should be applied   
    padFactHilbert = 1; % factor by which signal is padded before Hilbert trafo
    isEquidistGrid = true; % boolean determining if the image grid is equidistant    
    figHandle = 'last'; % handle to figure showing the image (other options: [](no fig),'new','last') 
    figTitle = ''; % title of figure (if .figHandle=[]: no effect)
    supressDrawnow = false; % boolean deciding of figure is immediately plotted, set true to add more curves to figure;    
end

methods
    
%% CONTRUCTOR: 

    function obj = visParamHandler(profile)
        %FILTPARAMS Construct an instance of this class
        %   Detailed explanation goes here
        
        if nargin<1
           profile = 'fast'; 
        end
        
        % set parameters according to profile:
        switch profile          
            case 'fast'
                % default values are used in this case!
            case 'accurate'
                obj.padFactHilbert = 2;                
            case 'memorySave'
                obj.padFactHilbert = 2;                
            otherwise
                error(['"%s" is not a valid profile! ',...
                 'Must be ''fast'',''accurate'' or ''memorySave'''], profile);
        end
        
    end

%% SETTERS & GETTERS:  


%% DISPLAY: 


end

end