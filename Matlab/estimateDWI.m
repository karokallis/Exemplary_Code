function [calcDWI, calcDWINormalized] = estimateDWI(Cmaps,TE, csPCa, urineMask)
%% Created by Karoline Kallis 05/2023
% Method to estimage calcualted DWI using different TE values
% for TE90, solve for Cmaps
% compare predicted DWI TEmin to observed TEmin1 and TEmin2

load('./data/params.mat');
bvals = [0;50;800;1500;3000];
%bvals = [0;500;1000;2000];
if csPCa == 0
    %params.T2s= [40, 40, 40, 40];
    params.T2s= [7.64, 25.37, 64.26, 15.27];
else
    params.T2s = [20.07, 47.92, 22.41, 14.53];
end
calcDWI = zeros(size(Cmaps,1),size(Cmaps,2),size(Cmaps,3),numel(bvals));
for bvali=1:numel(bvals)
    for id1 = 1:numel(params.ModelADCs)
        calcDWI(:,:,:,bvali) = calcDWI(:,:,:,bvali) + Cmaps(:,:,:,id1).*exp(-bvals(bvali).*params.ModelADCs(id1)) .* exp(-TE./params.T2s(id1));
        if(bvals(bvali)==0)
            tmp = calcDWI(:,:,:,bvali);
            normbladder = median(tmp(urineMask>0.5),'all');
        end
    end
end
calcDWINormalized = calcDWI./normbladder;
end

