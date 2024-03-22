%% CALCULATEFLUENCECOMPENSATIONMAP calculates a map for fluence ompensation
%
% DESCRIPTION: 
% Use the ultrasound image to identify different materials by thresholding
% and then calcualte an exponentially decaying fluence estimate in each
% material using the absorption coefficients of the materials at the
% specific wavelengths!
%
% INPUTS:
% - z_axis: array with z positions of the pixels of the fluence compensation map [m]
% - x_axis: array with z positions of the pixels of the fluence compensation map [m] 
% - absCoeffs: e.g. matrix of absorption coefficients of the different materials
%              (e.g. tissue and water) at the two wavelengths [1/cm].  
% - imgUS: ultrasound image defined on the same axes as the map that can be 
%          used for differentiating between the materials (water and tissue) [a.u]
%
% OUTPUT:
% - fluenceCompensationMap: (Nz-by-Nx matrix that compensates for the
%                            fluence loss when being multiplied to an image [1]
%
% AUHTOR: 000 ENTER YOUR TEMNAME 000

function fluenceCompensationMap = calculateFluenceCompensationMap(z_axis,x_axis, absCoeffs, imgUS)

    %  create water-tissue mask
    % (this can be a binary mask that is 1 in each pixel that represents 
    %  tissue and 0 for each pixel that represents water. You can use the US 
    %  image and create a mask from it by thresholding)
    tissueMask = (imgUS > -30); %Everything above -30 is true (1) and not water

    %  create Fluence map:
    % (use your mask to calcualte the fluence at each depth
    fluenceMap = zeros(length(z_axis),length(x_axis)); % Matrix with zeros with same dimension image

    for x = 1:length(x_axis)
%        fluenceWater = zeros(length(z_axis));
%        fluenceTissue = zeros(length(z_axis));

        for z = 1:20
            fluenceWater = exp(-z_axis(~tissueMask(z,x)) * absCoeffs(1));
            fluenceTissue = exp(-z_axis(tissueMask(z,x)) * absCoeffs(2));
            disp(fluenceWater)
            disp(fluenceTissue)


            if (size(fluenceWater) == [1,1])
                fluenceMap(z,x) = fluenceWater;
            else
                fluenceMap(z,x) = fluenceTissue;
            end

        end

    end
    
    imagesc(fluenceMap)
    
    %  convert fluence map to fluence compensation map:
    %  (which is later multiplied with the PA images)
    fluenceCompensationMap = 1./fluenceMap;

end