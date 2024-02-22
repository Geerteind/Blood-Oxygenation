%% APPLYPARECONSTRUCTION reconstructs photoacoustic receive data into image data
%
%USAGE:
%  imgData = applyPAReconstruction(ReceiveData, fs,c0, x_elem, z_axis,x_axis, FNumber)
%
%INPUTS:
%  - ReceiveData: (Nt-by-Nel-by-Nfr) array containing PA receive data [a.u.]
%  - fs: scalar sample frequency of imaging system / receive data [1/s]
%  - c0: scalar speed of sound used for the reconstruction [m/s]
%  - x_elem: vector with x-positions of transdcuer elements [m] (z-positions are assumed 0)
%  - z_axis: vector with z-positions of the image rows of the output imgData [m]
%  - z_axis: vector with x-positions of the image columns of the output imgData [m]
%  - FNumber: positive scalar that defines the ratio of pixel depth to active aperture [1]
%
%OUTPUTS:
%  - imgData: (Nz-by-Nx) array with reconstructed RF image [a.u.]
%
%AUTHOR: 000 FILL IN YOUR TEAMNAME 000

function imgData = applyPAReconstruction(ReceiveData, fs,c0, x_elem, z_axis,x_axis, FNumber)

    %  get other parameters (use these parameters to calculate the open boxes):
    dt  = 1/fs;                 % temporal spacing between samples [s]
    Nt  = size(ReceiveData,1);  % number of temporal samples 
    Nel = size(ReceiveData,2);  % number of transducer elements
    Nfr = size(ReceiveData,3);  % number of acquired frames

    %  initialize image:
    imgData = zeros(length(z_axis), length(x_axis), Nfr);

    tic;
    %  reconstruct all frames:
    for i_fr = 1:Nfr
        
        %  loop over all pixels:
        for i_z = 1 : 000 
            for i_x = 1 : 000 
                
                %  sum over all elements:
                for i_el = 1 : 000     

                    %  get pixel position and element position for this iteration:
                    z_pix = 000 ;
                    x_pix = 000 ;
                    x_el  = 000 ;

                    %  restrict F-number (only add elements under a certain angle to the pixel):
                    if abs(z_pix/(x_el-x_pix)) < FNumber, continue; end
                    
                    %  compute PA delay as travel distance converted to a time [s]:
                    %  (according to your findings in the Research phase)
                    tau_delays = 000 ;

                    %  convert delay from seconds to time index:
                    %  (this is still a real number, not an integer)
                    i_delays = 000 ;
                    
                    %  convert (real) index into integer index:
                    i_delaysInt = 000 ;

                    %  don't add data for indexes that are out of range:
                    if (i_delaysInt<1)||(i_delaysInt>Nt), continue; end
                    
                    %  retrieve Data at the respective delays:
                    %  (you can read out the data at the rounded integer index,
                    %   but for better image quality, you can also apply a 
                    %   weighted sum of values at multiple indexes around 
                    %   the real index, e.g. by usning linear interpolation) 
                    imgData_pix = 000 ;

                    %  add to pixel value by summing over all channels:
                    imgData(i_z, i_x, i_fr) = 000 ;
                end
            end
        end

% uncomment plot command for debugging:
%        % display reconstructed image :
%         imagesc(x*1e3, z*1e3, abs(hilbert(imgData(:,:,i_fr)))); axis image; colormap gray;
%         xlabel('x [mm]'); ylabel('z [mm]'); drawnow ;

        toc;
    end

end