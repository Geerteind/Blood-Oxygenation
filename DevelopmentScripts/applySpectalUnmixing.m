%% APPLYSPECTRALUNMIXING converts wavelengths images to component images
%
%USAGE:
%  [componentData1,componentData2] = applySpectalUnmixing(imgData1,imgData2, absCoeff_component1,absCoeff_component2)
%
%INPUTS:
%  - imgData1: (Nz-by-Nx) array containing PA image at 1st wavelength (750nm) [a.u]
%  - imgData2: (Nz-by-Nx) array containing PA image at 2nd wavelength (850nm) [a.u]
%  - mu_component1: 2-element column array with absorption coefficents of 1st component (HB) [a.u.]
%  - mu_component2: 2-element column array with absorption coefficents of 2nd component (HBO2) [a.u.]
%
%OUTPUTS:
%  - componentData1: (Nz-by-Nx) array containing image of 1st component (HB) [a.u]
%  - componentData2: (Nz-by-Nx) array containing image of 2nd component (HBO2) [a.u]
%
%AUTHOR: 000 FILL IN YOUR TEAMNAME 000 

function [componentData1,componentData2] = applySpectalUnmixing(imgData1,imgData2, mu_component1,mu_component2)

    %  get dimensions of image:
    [Nz, Nx] = size(imgData1);

    %  reshape imgData into a Nwl-by-Nx*Nz array with all pixels of the
    %  first wavelength in the first row and all pixels of the second 
    %  wavelenght in the second row:
    imgData1flat = reshape(imgData1,1,Nz*Nx);
    imgData2flat = reshape(imgData2,1,Nz*Nx);
    spectralData = [imgData1flat,imgData2flat];
    
    %  create a mixing matrix:
    %  (The matrix consists of absorption coefficients and maps a
    %  column vector of component data [HB;HBO2] to a column vector of 
    %  wavelength data [750;850])
    mixMatrix  = [mu_component1; mu_component2];
    disp(mixMatrix)
    %  apply unmixing:
    %  (find an inverse solution to the mixing matrix. You can also look for  
    %   more advanced inversion approaches in literature! For example, you  
    %   can get a least-squares solution by using the backslash operator "\", 
    %   or a solution that enforces positive values using "lsqnoneg") :
%    componentData = mixMatrix \ spectralData;
    componentData = zeros(2,Nz*Nx) ;
    
    for i=1:Nx*Nz
        componentData(:,i) = lsqnonneg(mixMatrix, spectralData(:,i));
    end
    
    %  reshape the (blood) component data from a two-row matrix with all 
    %  pixels of one component per row into a 3D array of two image matrixes
    %  with the two components in the third dimension:
    componentData1 = reshape(componentData(1,:), Nz, Nx);     
    componentData2 = reshape(componentData(2,:), Nz, Nx);  
    
end

