% loadDataSimulation loads simulation receive data and metadata from a mat file
%
% USAGE:
% [RXD750,RXD850,fs,c0,x_elem] = loadDataSimulation(filePath)
%
% INPUTS:
% - filePath: string with relative or absolute path to the matfile to read
%
% OUTPUTS:
% - RXD750: 2D array (t,x) containing PA receive data at 750nm [a.u]
% - RXD850: 2D array (t,x) containing PA receive data at 750nm [a.u]
% - fs:     scalar sample frequency of receive dataset [1/s]
% - c0:     scalar speed of sound used for simulation [m/s]
% - x_elem: vector with x-positions of simulated transducer elements [m]
%
% AUTHOR: h.schwab@tue.nl (Feb22)

function [RXD750,RXD850,fs,c0,x_elem] = loadDataSimulation(filePath)

    % settings:
    upsampleFactor = 4;

    % load mat file:
    data = load(filePath);

    % get receive data:
    RXD750 = interpft(data.channelData{1}, size(data.channelData{1},1)*upsampleFactor, 1);
    RXD850 = interpft(data.channelData{2}, size(data.channelData{2},1)*upsampleFactor, 1);

    % get metadata:
    fs     = data.sample_frequency * upsampleFactor;
    x_elem = data.Elem_position;% x_elem = data.Elem_position(1:2:end-1);% (-0.02016 : 3.1500e-04 : 0.019845);
    c0     = data.sos;

end