%RECONMEX applies a fast 2D reconstruction using mex code
%
%INPUT:
% - RXData: 2D or 3D array with receive data with temporal samples
%           of all elements in columns and different transmission
%           events in the third dimension, which will be compounded.
% - Z: 1D array or 2D grid with z-positions of the volume pixels [m]
% - X: 1D array or 2D grid with x-positions of the volume pixels [m]
% - dt: temporal spacing of RXData [s]. If .filter was applied prior to this
%       method, the temporal spacing might not equal 1/fs, because the 
%       data might have been upsampled to increase the acuracy of DAS. 
% - recompileMex: boolean deciding if the mex file should be
%                 recompiled before execution. The file will be
%                 deployed in the mexfiles directory.
%OUTPUT:
% - imgData: compounded high-frequency image dataset on the grid defined 
%            by z and x (no envelope detection, no log-compression). Can be
%            RF or IQ, depending on if "RXData"  is RF or IQ.
%NOTE:
% The method used the volume reconstruction mexcode "recon3Dmex" and
% only evaluates it in the y=0 plane. For better performance, the
% mexcode shouldbe rewritten for 2D reconstruction!
%
%AUTHOR: hans-martin.schwab@web.de (Jan22) (Mexcode by JW Muller)

%TODO:
% - write 2D mex code
% - implement apodization!

function imgData = reconMex(obj, RXData_filt, dtu, z, x, recompileMex)
    if nargin<=1, compileMex(); return; else, narginchk(5,6); end
    t_rec = tic;       
    % input checking:
    if nargin<6, recompileMex = false; end 
    obj.checkConsistency(RXData_filt); 
    %can't find check for gpuArray??? So just gather it        
    RXData_filt = gather(RXData_filt);
    x = gather(x);
    y = 0;
    z = gather(z);
    dtu = gather(dtu);
    % check if input axes are vectors of grids and turn into grids:
    if any(size(x) == 1), [z,x,y] = ndgrid(z,x,y); end         
    % assure input data is complex even for RF:
    if isreal(RXData_filt), isReal = true; RXData_filt = complex(RXData_filt); 
    else, isReal = false; end 
    % define 3D RX aperture:
    elemPos3D = [obj.elemPos,zeros(obj.Na,1)];
    % consider optional TX aperture if defined:
    if isempty(obj.elemPosTX), elemPositionTX3D = elemPos3D;
    else, elemPositionTX3D = [obj.elemPosTX,zeros(obj.Na,1)]; end
    % set interpolaoitn method:
    switch obj.reconParams.interpMethod
        case 'next', interpMethod = 1;
        case 'linear', interpMethod = 2;
        case 'phaseLinear', interpMethod = 3;
        otherwise, error('This .reconParams.interpMethod is not supported (must be on of: "next","linear","phaseLinear")!');
    end
    % set additional input arguments:
    offset = obj.TXPulseDelays;  
    % compile mex file if necessary:
    if recompileMex, compileMex(); end
    % apply reconstruction:
    imgData =  recon3Dmex(RXData_filt,...
                          elemPos3D,...
                          z,...
                          x,...
                          y,...
                          obj.c0,...
                          obj.fc,...
                          1/dtu,...
                          obj.TXPulseDelays,... 
                          offset,...
                          obj.isPA,...
                          obj.reconParams.applyCompounding,...
                          interpMethod,...
                          elemPositionTX3D);
    if isReal, imgData = real(imgData); end
    % store computation time:
    obj.compTimeRecon = toc(t_rec);
end

    %% Local function:
    
function compileMex()
    % go to directory where mex reconstruction is stored 
    prevPath = pwd; 
    curPath = mfilename('fullpath'); 
    curPath = [curPath(1:end-(18)),'\mexfiles'];  
    cd(curPath); % mex file is stored in same directory as Recon3D 
    mex -R2018a  recon3Dmex.cpp CXXFLAGS="$CXXFLAGS -fopenmp -march=native" CXXOPTIMFLAGS="-Ofast -DNDEBUG" LDFLAGS="$LDFLAGS -fopenmp";
    cd(prevPath); %go back to previous path 
end
    