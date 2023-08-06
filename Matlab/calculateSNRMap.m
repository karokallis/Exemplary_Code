function [SNRmap , OutlierMap]= calculateSNRMap(data, prostate)
% Created by Karoline Kallis 6/1/2023
% function to calculate the SNR map for the prostate, outside of prostate
% considered noise

% input: data =  5D matrix - x-y-z-cmaps-patients
%        prostate = 5D matrix - x-y-z-cmaps-patients -- repeated
%        segmentation for each cmap

% output: SNRmap for each compartment and each patient
         %x-y-z-cmaps-patients
         
SNRmap = zeros(size(data));
OutlierMap = zeros(size(data));

for pati=1:size(data,5)
    for cvali=1:size(data,4)
        tmpdata = squeeze(data(:,:,:,cvali,pati));
        tmpprostate = squeeze(prostate(:,:,:,cvali,pati));
        
        % calcualte Noise
        % Noise is RMSE of outside area of prostate
        
        noise = sqrt(mean((tmpdata(tmpprostate>0.5)-mean(tmpdata(tmpprostate>0.5))).^2));
        SNRmap(:,:,:,cvali,pati) = tmpdata./noise;
        
        prc = prctile(tmpdata(tmpprostate<0.5),99);
        OutlierMap(:,:,:,cvali,pati) = tmpdata./prc;
        
    end
end

end

