function C_vol = calculateRSI_TE(DWIdata,ubvals, TE, csPCa)
%Karoline Kallis 4/18/2023 function calculate RSI including TE and T2
%effect
% INPUT: 
% DWIdata = noise corrected and avareged DWI image
% ubvales = unique b-values acquired
% TE for this particular acquistion


% Sobs 1*n-bvalues array of signal intensity value -- N*M*3 DWI transformed
% to 1D where each row represents a different b-value acquistion

% % Compartment 1
% estT2_WP = 7.64;
% estT2_cancer = 20.07;
% 
% % Compartment 2
% estT2_WP = 25.37;
% estT2_cancer = 47.92;
% 
% % Compartment 3
% estT2_WP = 64.26;
% estT2_cancer = 22.41;
% 
% % Compartment 4
% estT2_WP = 15.27;
% estT2_cancer = 14.52;
% 
load('./data/params.mat');
if csPCa == 0
    params.T2 = [7.64, 25.37, 64.26, 15.27];
else
    params.T2 = [20.07, 47.92, 22.41, 14.53];
end

x(1,:) = params.ModelADCs;
x(2,:) = params.T2;

Sobs = zeros( length(ubvals), numel(DWIdata(:,:,:,1)) );
for b = 1:length(ubvals)
    vol_b = DWIdata(:,:,:,b);
    Sobs(b,:) = vol_b(:)';
end
temp_data = struct('Sobs', Sobs, 'bvals', ubvals', 'TE', TE);

[~, aic, C_mat] = fit_RSI_model_mTE(x, temp_data);

C_vol = zeros(size(DWIdata,1), size(DWIdata,2), size(DWIdata,3), size(C_mat,1));
for c = 1:size(C_mat, 1)
    C_vol(:,:,:,c) = reshape( C_mat(c,:), size(C_vol,1), size(C_vol,2), size(C_vol,3) );
end





end

