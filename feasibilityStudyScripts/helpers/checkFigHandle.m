% function [] = checkFigHandle(figHandle)
%
% checks the validity of a figure. If the figure handle exists, but the
% figure has been closed, a warning is printed and a new figure is created.
% Also, if the figHandle is not a figure handle, a new figure is created. 
% Otherwise the figure belonging to the handle becomes the current figure.
%
% author: hans-martin.schwa@web.de (May19)

function isValidFig = checkFigHandle(figHandle)

    % initilaize as false
    isValidFig = false;
    
    if ~isempty(figHandle)      
        % set to true (plot a figur ein any case)
        isValidFig = true;        
        % create new figure if figHandle is no handle:
        if ~ishandle(figHandle)
            %create new figure if figHandle matches string 'new':
            % (otherwise no new figure is created and hence the last opened
            %  figure is used for the plot)
            if strcmp(figHandle,'new')
                figure;
            end
        else               
            % create new figure if figHandle is not valid (was closed):        
            if ~isvalid(figHandle)
                warning(['params.figureHandle is not a valid figure handle.',...
                         'The figure might have been closed. A new figure is created.']);
                figure;   
            else
                % open figHandle related figure:
                figure(figHandle);                 
            end           
        end
    end    
    
end