function [pred, aic, beta] = fit_RSI_model_mTE(x, temp)
% created by Anders and adpated by Chris 
% temp is a structure with the fields:
% x(1,:)= ADC values for each compartment
% x(2,:) = T2 value for each compartment

% bvals - unique b-values
% TE - Echo time
% Sobs - flattened signel intensity values for each b-value

% returns: pred = prediction DWI (similar to synthetic b-values)
% aic = defines the goodness of the fit
% beta = cmaps -- called beta because of brain people stuff


if ~isfield(temp, 'bvals')
    error('%s - Error! bvals is not defined!\n', mfilename);
end
if ~isfield(temp, 'TE')
    error('%s - Error! TEs is not defined!\n', mfilename);
end
if ~isfield(temp, 'Sobs')
    error('%s - Error! Sobs is not defined!\n', mfilename);
end
% if ~isfield(temp, 'tol')
%     temp.tol = 0;
% end

bvals = temp.bvals; 
TE = temp.TE;
Sobs = temp.Sobs; 
% tol = temp.tol;

ADCs = x(1,:);
T2s = x(2,:);

A = [];
%B = [];
for id1 = 1:numel(ADCs)
    %B = cat( 2, B, reshape( exp(-bvals.*ADCs(id1)), [], 1));
    A = cat( 2, A, reshape( exp(-bvals.*ADCs(id1)) .* exp(-TE./T2s(id1)), [], 1));
end

% if 0 % Old, slow version
%     AtA = transpose(A)*A;
%     Ainv = inv(AtA+tol*mean(diag(AtA))*eye(size(AtA)))*transpose(A);
%     beta = Ainv*Sobs;
%     negind = find(any(beta<0,1));
%     negind = find(any(true(size(beta)))); % Do lsqnonneg for all -- check if different from standard pseudoinverse for any with non-negative values
%     tic
%     for id1 = 1:numel(negind)
%         b_temp = Sobs(:,negind(id1));
%         beta(:,negind(id1)) = lsqnonneg(A, b_temp, optimset('Display', 'off'));
%     end
%     toc;
%     beta_bak = beta;
% else
    beta = lsqnonneg_amd(A,Sobs);
% end

pred = A*beta;
%cost = sum(mean((pred - Sobs).^2, 1)./mean(Sobs,1), 2);
rss = sum((pred - temp.Sobs).^2, 1);
numComp = length(x);
aic = 2*numComp + 4*log(rss./4) + (2*numComp*(numComp+1))./(4-numComp-1); % Not sure this is correct
% 
% if 0
%     resmat_bak = A*beta_bak - Sobs;
%     resmat = A*beta - Sobs;
%     figure(666); subplot(2,2,1); plot(beta(1,:)-beta_bak(1,:)); subplot(2,2,2); plot(beta(2,:)-beta_bak(2,:)); subplot(2,2,3); plot(beta(3,:)-beta_bak(3,:)); subplot(2,2,4); plot(sqrt(sum(resmat.^2,1))-sqrt(sum(resmat_bak.^2,1)));
%     [sv si] = sort(sqrt(sum(resmat.^2,1))-sqrt(sum(resmat_bak.^2,1)),'descend');
%     ii = si(1);
%     beta(:,ii)'
%     beta_bak(:,ii)'
%     beta_tmp = A(:,1:2)\Sobs(:,ii); % Works -- change code to fit for all combinations of non-zero
% end

end
