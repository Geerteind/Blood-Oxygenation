%RECONPARAMHANDLER Class to manage reconstruction parameters of DASReconstructor class
%
%PROPERTIES:  
%.apodCutoffAngle    = 20;         maximum supported angle between element normal and pixel [deg] within the interval [0,90]
%.apodWeights        = [];         can be an array with weights between cutoff or [] for binary apodization or 'hanning' 
%.interpMethod       = 'next';     string with interpoalation method used in DAS {'next','linear','phaseNext','phaseLinear'}
%.insoniWeightType   = 'WLBinary'; type of weighting the insonified area (options:'WLBinary','WLConti','off')
%.insoniWeightParams = .1;         parameters for insonification weigths 
%.applySerialFrameComp  = false;   applies reconstruction in serial over frames to avoid GPU timeout issues
%.applyCompounding   = true;       If true, the result of .reconstruct is compounded, if false the third dimension equals .Nacq
%.acqWeights         = [];         applies weights to the TX events in RXData, e.e. to apply synthetic TX apodization
%.figHandle          = [];         handle to figure showing the elements and the pixels (if []: no figure) 
%
%AUTHOR: 
% Hans-Martin Schwab (hans-martin.schwab@web.de)

classdef reconParamHandler < matlab.mixin.Copyable
    %RECONPARAMHANDLER Class to manage reconstruction parameters of DASReconstructor class
    
properties
    
    % maximum supported angle between element normal and pixel [deg] within the interval [0,90]
    apodCutoffAngle = 20;
    % string with interpoalation method used in DAS {'next','linear','phaseNext','phaseLinear'}
    interpMethod = 'next';     
    % can be an array with weights between cutoff or [] for binary apodization or 'hanning'
    apodWeights = []; 
    % type of weighting the insonified area (options:'binary','continous')
    insoniWeightType {mustBeMember(insoniWeightType,{'off','WLBinary','WLConti'})}  = 'off'; 
    % parameters for insonification weigths
    insoniWeightParams = .1;
    % applies reconstruction in serial over frames to avoid GPU timeout issues:
    applySerialFrameComp  = false;
    % applies compounding of all acquisitions:
    applyCompounding = true;
    % applies weights to the TX events in RXData, e.e. to apply synthetic TX apodization:
    acqWeights = [];
    % handle to figure showing the elements and the pixels (if []: no figure)
    figHandle = [];
        
    
    % uncontinued and only still in the property list to trwo an error,
    % because the property was moved form reconParams to the Recon object
    % itself:
    normalAngleRX = [];
    
end

methods
    
%% CONTRUCTOR: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    function obj = reconParamHandler(profile)
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
                obj.apodCutoffAngle = 30;
                obj.apodWeights = 'hanning';
                obj.insoniWeightType = 'WLConti';
                obj.insoniWeightParams = 1; 
            case 'memorySave'
                obj.applySerialFrameComp = true;               
            otherwise
                error(['"%s" is not a valid profile! ',...
                 'Must be ''fast'',''accurate'' or ''memorySave'''], profile);
        end
        
    end
    
%% HELPERS:

%     % retrieve cutoff angle by detmeirning attenuation of sensitivity curve 
%     function setMaxSensiAtten(obj,atten, elemWidth, fc)
%         obj.
%     end

%     % retrieve cutoff angle by detmeirning attenuation of sensitivity curve 
%     function setFNumber(obj, FNumber)
%         obj.
%     end

%% SETTERS & GETTERS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    function set.apodCutoffAngle(obj, val)
        if val == 90, val = 89.999; end
        assert(val<90, 'An cutoff angle higher than 90 degrees cannot be set!')
        assert(val>0, 'The cutoff angle must be >=0!');
        obj.apodCutoffAngle = val;
    end

    function set.interpMethod(obj, val)
        % list of valid interpolation methods:
        validInterpMethods = {'next','linear','phaseNext','phaseLinear'};
        % interp method can be provided as number and is then converted to
        % a string:
        if isnumeric(val)
            assert(round(val)>0&&round(val)<=length(validInterpMethods),...
                   [ 'Cannot evaluate reconParam: .interpMethod ! '...
                    ,' If method is provided as numeric value it must be between 1 and '...
                    ,num2str(length(validInterpMethods)),'!'])
            obj.interpMethod = validInterpMethods{val};
        % or it is already a string:
        else
            assert(any(strcmp(val,validInterpMethods)), [ 'Cannot evaluate reconParam: .interpMethod! '...
                   ,' Must be one of the options: {"next","linear","phaseNext","phaseLinear"}!'])
            obj.interpMethod = val;
        end
    end    
    
    % uncontinued proptery:
    function set.normalAngleRX(obj, val)
       error(['The property ".normalAngleRX" was moved from the ',...
              '"ReconParamHandler" class to the "Recon2D" / "Recon3D" class. ',...
              'Please change all calls of "obj.reonParams.normalAngleRX" to ',...
              '"obj.normalAngleRX"!']);
    end
    function val = get.normalAngleRX(obj)
       error(['The property ".normqlAngleRX" was moved from the ',...
              '"ReconParamHandler" class to the "Recon2D" / "Recon3D" class. ',...
              'PLease change all calls of "obj.reonParams.normalAngleRX" to ',...
              '"obj.normalAngleRX"!']);
    end
    
%% DISPLAY: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

end
