function [receiveDataRF, fs, x_elem] = loadDataAcousticX(filePath)


    fid = fopen(filePath, 'rb');
    
    try 
        data = fread(fid, 'int16');
    catch ME
        fclose(fid);
        rethrow(ME);
    end
    
    frame_num = size(data,1)/1024/128;
    receiveDataRF = reshape(data, [1024 128 frame_num]);
    
    fs = 60e6;
    x_elem = ((0:127)-64) * 2.2e-4;

end