function [receiveDataRF, fs, x_elem] = loadDataAcousticX(filePath, dataType)
    
    % input check:
    assert(any(strcmp(dataType,{'US','PA'})), 'Input argument "dataType" must be "US" or "PA"!');

    % settings:
    if strcmp(dataType,'PA')
        upsampleFactor = 4; % positive integer, by which dataset is upsampled
        fs = 40e6; % sample frequency of loaded data [Hz]
    else
        upsampleFactor = 1; % positive integer, by which dataset is upsampled 
        fs = 20e6; % sample frequency of loaded data [Hz]
    end
    
    % set further system property outputs:
    x_elem = ((0:127)-64) * 0.315e-3; % element positions in x-direction (along aperture)
    
    % load samples from binary file:
    fid = fopen(filePath, 'rb');
    try 
        data = fread(fid, 'int16');
    catch ME
        fclose(fid);
        rethrow(ME);
    end
    
    % retrieve RF data:
    frame_num = size(data,1)/1024/128;
    receiveDataRF = flip(reshape(data, [1024 128 frame_num]),2);
       
    % some extra processing for PA data:
    if strcmp(dataType,'PA')
        
        % correct for laser delay:
        receiveDataRF(1:25,:,:) = [];
        
        % reduce electrical distorption by removing average of first 200 samples
        receiveDataRF(1:200,:,:) = receiveDataRF(1:200,:,:) - mean(receiveDataRF(1:200,:,:),3);
        % receiveDataRF(1:200,:,:) = 0;
        
    end
    
    % upsample dataset if upsample factor is larger than 1:
    if upsampleFactor>1
        
        % adapt sample frequency:
        fs = fs * upsampleFactor;
        
        % pad zeros to data, upsample, and remove padding:
        receiveDataRF = cat(1, receiveDataRF, zeros(size(receiveDataRF)));
        receiveDataRF = interpft(receiveDataRF, size(receiveDataRF,1)*upsampleFactor, 1);
        receiveDataRF = receiveDataRF(1:end/2,:,:);
        
    end

end