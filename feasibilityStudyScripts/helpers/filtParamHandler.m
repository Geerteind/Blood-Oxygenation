%FILTPARAMHANDLER Class to manage filter parameters of DASReconstructor class
%
%PROPERTIES:
%.upsampleFactor = 1;          % factor by which signal is upsampled for better interpolation
%.type           = 'simpleBP'; % string (cell array?) with param type (options:'off','antialias','simpleBP','manualTD','manualFD') 
%.params         = [1,2];      % actual parameters of filters (cell array?)
%.inputFormat    = 'auto';     % string defining format into which input is converted (options:'auto','RF','IQ','IQDemod')
%.outputFormat   = 'IQ';       % string determining desired output format (options:'auto','RF','IQ','IQDemod')  
%.figHandle      = [];         % handle to figure showing the spectra of data and filter (if []: no figure)
%.attenCoeff     = 0;          % scalar determining the attenuation coefficient of the center frequency in dB/cm used for compensation
%AUTHOR: 
% Hans-Martin Schwab (hans-martin.schwab@web.de)

classdef filtParamHandler < matlab.mixin.Copyable

properties (SetAccess = private)
   
    validTypes    = {'off','antialias','simpleBP','manualTD','manualFD'}; % cell array of strings of supported filter types
    defaultParams = struct('off', [],...
                           'antialias', [1,1,20],...
                           'simpleBP', 1,...
                           'manualTD', fir1(40,[.15,.85]),...
                           'manualFD', hanning(100)); 
end    
    
properties
    
    % factor by which signal is upsampled for better interpolation
    upsampleFactor = 4;
    % string (cell array?) with param type (options:'off','antialias','simpleBP','fir','manual') 
    filtType = 'simpleBP';
    % actual parameters of filters (cell array?)
    params = 1;
    % string defining format into which input is converted (options:'auto','RF','IQ','IQDemod')
    inputFormat = 'auto';
    % string determining desired output format (options:'auto','RF','IQ','IQDemod')  
    outputFormat = 'IQ'; 
    % handle to figure showing the spectra of data and filter (if []: no figure)
    figHandle = []; 
    % scalar determining the attenuation coefficient of the center frequency in dB/cm used for compensation:
    attenCoeff = 0;
    
end

methods
    
%% CONTRUCTOR: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    function obj = filtParamHandler(profile)
        %FILTPARAMHANDLER Class to manage filter parameters of DASReconstructor class
        
        if nargin<1
           profile = 'fast'; 
        end
        
        % set parameters according to profile:
        switch profile          
            case 'fast'
                % default values are used in this case!
            case 'accurate'
                obj.upsampleFactor = 8;
            case 'memorySave'
                
            otherwise
                error(['"%s" is not a valid profile! ',...
                 'Must be ''fast'',''accurate'' or ''memorySave'''], profile);
        end
    end

%% SETTERS & GETTERS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % type
    function set.filtType(obj, val)
       assert(ischar(val),'property "type" must be a string!');
       assert(ismember(val,obj.validTypes),'filter type is not supported!');
       obj.filtType = val;
%        obj.resetParams;
    end
    
    % params
    function set.params(obj, val)
        
       if ismember(obj.filtType,{'simpleBP','antialiasing','manualTD'}) %QUICK AND DIRTY SOLUTION BACUSE 'manualFD' can be funhandle        
       assert(strcmp(class(val),class(obj.defaultParams.(obj.filtType))),...
              'filter parameter format is not correct for this filter type!');
       end
       
       % check number of elements for each filter individually
       % TODO
       
       % apply parameters (CHECK EACH VALUE FOR EACH TYPE!)
       obj.params = val;
    end
    
%% HELPERS:

%NOT COMPATIBLE WITH CURRENT IMPLEMENTATION:
%     % reset to default filter params:
%     function resetParams(obj)
%         % set params to default params of current filter type:
%         obj.params = obj.defaultParams.(obj.filtType);
%     end
    
%% DISPLAY: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

end

