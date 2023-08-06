function [fout,fout_intercept, R2out,f,dataTilde] = estimateScalingFactor(reference,data, repetition, batch, plotFigures)
%Created by Karoline Kallis 05/25/2023
% Function estimate scaling factor between two datasets with different TE
% using linear regression
% Idea: is to calcualte linear scaling factor to transform TE90 series to
% TEmin acquistion

% input = Cmaps normalized to median urine of B0 *100
% reference = cell of all datasets with all Cmaps
% reference.data.outmaps_normalized

% data = cell of all datasets with all Cmaps of different TE time e.g.,
% TE90

% output fout = hopefully a linear scaling factor between the datasets
% f = all output scalers for repeated regression system
% dataTilde = f*data supposedly calibrated data


% 1) Create matrix for all patients
% Under the assumption that all outmaps are sampled in the same manner!
% Potentital need to add a resampling pitfall for this case
dataM = zeros([size(data.outmaps{1}) length(data.outmaps)]);
referenceM = zeros([size(reference.outmaps{1}) length(reference.outmaps)]);
%SNRmap = zeros([size(reference.outmaps{1}) length(reference.outmaps)]);
dataTilde = zeros([size(data.outmaps{1}) length(data.outmaps)]);
if(size(dataM)~=size(referenceM))
    disp('Error: Matrices have a different size!!!')
end

ProstateSegmentation = zeros([size(reference.outmaps{1}) length(reference.outseg)]);
f = zeros(size(referenceM,4),repetition);
R2 = zeros(size(referenceM,4),repetition);


    for pati=1:length(data.outmaps)
        dataM(:,:,:,:,pati) = data.outmaps_normalized{pati};
        referenceM(:,:,:,:,pati) = reference.outmaps_normalized{pati};
        %     pati
        for cvali=1:size(referenceM,4)
            ProstateSegmentation(:,:,:,cvali,pati) = reference.outseg{pati};
        end
        
        if((sum(dataM(:,:,:,:,pati),'all')==0 || sum(referenceM(:,:,:,:,pati),'all')==0 || sum(ProstateSegmentation(:,:,:,:,pati),'all')==0))
            disp('Dataset empty for patient: ', pati)
        end
        
    end
    
    %% eliminate voxel which are too noisy
    [SNRmapData, ~]= calculateSNRMap(dataM, ProstateSegmentation);
    [SNRmapReference, ~] = calculateSNRMap(referenceM, ProstateSegmentation);
    
    ProstateSegmentation(SNRmapData<1.0|SNRmapReference <1.0) = 0;
    % ProstateSegmentation(OutlierData>1.0|OutlierReference >1.0) = 0;
for repi=1:repetition 
%     ProstateSegmentation = origSeg;
    patRadi = randperm(length(data.outmaps),min(batch, length(data.outmaps)));
%     for leavi = patRadi
%         ProstateSegmentation(:,:,:,:,leavi)=0;
%     end
    %% Test if linear regression is adequad
    % estimate each function f(x) that each feature has a linear relationship
    % to the target
    % y = X*beta
    % Temin  = slope*Te90 + intercept
    
    % 1.) linearity
    
    for cvali = 1:size(referenceM,4)
        %     tmpData = squeeze(dataM(:,:,:,cvali,:).*ProstateSegmentation(:,:,:,cvali,:));
        %     tmpData = sqrt(squeeze(max(tmpData,[],[1,2,3])).^2);
        %
        %     tmpReference = squeeze(referenceM(:,:,:,cvali,:).*ProstateSegmentation(:,:,:,cvali,:));
        %     tmpReference = sqrt(squeeze(max(tmpReference,[],[1,2,3])).^2);
        %     tmpReference = sqrt(mean(squeeze(referenceM(:,:,:,cvali,:),[1,2,3])).^2);
        %     tmpProstate  = ProstateSegmentation(:,:,:,cvali,:);
        %     dataM(SNRmapData<1) = NaN;
        %     referenceM(SNRmapReference<1) = NaN;
%                     tmpData = sqrt(mean(squeeze(dataM(:,:,:,cvali,:)),4).^2);
%                     tmpReference = sqrt(mean(squeeze(referenceM(:,:,:,cvali,:)),4).^2);
        
        
        %     tmpData1 = sqrt(prctile(squeeze(dataM(:,:,:,cvali,:)),75,4).^2);
        %     tmpReference1 = sqrt(prctile(squeeze(referenceM(:,:,:,cvali,:)),75,4).^2);
        %
        %         tmpData2 = sqrt(prctile(squeeze(dataM(:,:,:,cvali,:)),95,4).^2);
        %         tmpReference2 = sqrt(prctile(squeeze(referenceM(:,:,:,cvali,:)),95,4).^2);
        
        tmpData = dataM(:,:,:,cvali,:);
        tmpReference = referenceM(:,:,:,cvali,:);
        %     %
        % I need to think of something better -- That is probably not
        % accurate at all
%         tmpProstate  = sum(squeeze(ProstateSegmentation(:,:,:,cvali,:)),4);
        tmpProstate = ProstateSegmentation(:,:,:,cvali,:);
        if ~isempty(patRadi)
            for leavi = patRadi
                tmpProstate(:,:,:,:,leavi) = 0;
            end
        end
        %     median(tmpData(tmpProstate>0.5),'all')
        %     max(tmpData(tmpProstate>0.5),[],'all')
        %
        %     median(tmpReference(tmpProstate>0.5),'all')
        
        if(plotFigures)
            x  = max(tmpReference(tmpProstate>0.5),[],'all');
            hv_C = linspace(0.001,x,1001);
            
            
            
            % Create cumulative histogram
            %       [hcD ~] = histcounts(tmpData,hv_C); % hv_C random binning // hc are counts // make histogram of voxel
            
            [hcD ~] = histcounts(tmpData(tmpProstate>0.5),hv_C); % hv_C random binning // hc are counts // make histogram of voxel
            %     figure
            %     histogram(tmpData(tmpProstate>0.5))
            hcD = hcD/sum(hcD);
            chcD = [0,cumsum(hcD)];
            %     [hcR ~] = histcounts(tmpReference,hv_C); % hv_C random binning // hc are counts // make histogram of voxel
            
            [hcR ~] = histcounts(tmpReference(tmpProstate>0.5),hv_C); % hv_C random binning // hc are counts // make histogram of voxel
            %     figure
            %     histogram(tmpReference(tmpProstate>0.5))
            hcR = hcR/sum(hcR);
            chcR = [0,cumsum(hcR)];
            
            figure;
            plot(hv_C, [0,hcD],'.r')
            hold on
            plot(hv_C, [0,hcR],'.b')
            title(['Histogram c' num2str(cvali) '-map'])
            legend('Data', 'Reference')
            hold off
            
            figure;
            plot(hv_C,chcD,'.r')
            hold on
            plot(hv_C,chcR,'.b')
            legend('Data', 'Reference')

            ylabel('cdf')
            title(['cdf c' num2str(cvali) '-map'])
            hold off
            
            %     figure
            %     plot(hv_C,log(chcD),'.r')
            %     hold on
            %     plot(hv_C,log(chcR),'.b')
            %     ylabel('log(cdf)')
            %     title(['cdf c' num2str(cvali) '-map'])
            %     hold off
            %
            %     figure
            %     plot(hv_C,1-log(chcD),'.r')
            %     hold on
            %     plot(hv_C,1-log(chcR),'.b')
            %     ylabel('1- log(cdf)')
            %     title(['cdf c' num2str(cvali) '-map'])
            %     hold off
            
            %     figure
            %     plot(tmpData(tmpProstate>0.5),tmpReference(tmpProstate>0.5),'.k')
            % %      plot(tmpData,tmpReference,'.k')
            % %     plot(1-log(tmpData(tmpProstate>0.5 & tmpReference>0.001)),1-log(tmpReference(tmpProstate>0.5 & tmpReference>0.001)),'.k')
            %     title(['c' num2str(cvali) '-map'])
            %     xlabel('log(data)')
            %     ylabel('log(ref)')
            
            fig = figure;
            qqplot(tmpData(tmpProstate>0.5),tmpReference(tmpProstate>0.5))
            %     qqplot(tmpData,tmpReference)
            %     qqplot(1-log(tmpData(tmpProstate>0.5 & tmpReference>0.001)),1-log(tmpReference(tmpProstate>0.5 & tmpReference>0.001)))
            title(['QQPlot c' num2str(cvali) '-map'])
            print(fig, ['images/QQplot_TEmin_TE90_SNR_batch',num2str(cvali)],'-dpng','-r100');
%             mdlCross = fitrlinear(tmpData(tmpProstate>0.5),tmpReference(tmpProstate>0.5),'CrossVal','on')
            %mdl1 = fitlm(rmoutliers(tmpData(tmpProstate>0.5),'mean'),rmoutliers(tmpReference(tmpProstate>0.5),'mean'))
        end
        mdl = fitlm(tmpData(tmpProstate>0.5),tmpReference(tmpProstate>0.5))
        %[breg,~,~,~,stats] = regress(tmpReference(tmpProstate>0.5),[tmpData(tmpProstate>0.5),ones(length(tmpData(tmpProstate>0.5)),1)])

        [breg,~,~,~,stats] = regress(tmpReference(tmpProstate>0.5),tmpData(tmpProstate>0.5));
        %b = tmpData(tmpProstate>0.5)\tmpReference(tmpProstate>0.5)
        %     mdl = fitlm(tmpData,tmpReference)
        %     mdl = fitlm(1-log(tmpData(tmpProstate>0.5 & tmpReference>0.001)),1-log(tmpReference(tmpProstate>0.5 & tmpReference>0.001)))
        %b = pinv(log(tmpData(tmpProstate>0.5 & tmpReference>0.001)))*log(tmpReference(tmpProstate>0.5 & tmpReference>0.001));
        %mdl.intercept
        
        if(plotFigures)
            fig = figure;
            plot(mdl)
            hold on
            plot(tmpData(tmpProstate>0.5),tmpData(tmpProstate>0.5)*breg(1),':m')
            title(['c' num2str(cvali) '-map'])
            %print(fig, ['images/Model_TEmin_TE90_b1_SNR_batch',num2str(cvali)],'-dpng','-r100');
        end
        % %     x = [0:1:100];
        % %     y1 = mdl.Coefficients.Estimate(2)*x+mdl.Coefficients.Estimate(2);
        % %     y2 = breg*x;
        % %     y3 = b*x;
        
        %     figure
        %     plot(x,y1,':r')
        %     hold on
        %     plot(x,y2,'-b')
        %     plot(x,y3,'.m')
        %     hold off
        
        % %% RANSAC
        % % RANSAC is accomplished with the following steps
        % %
        % % Randomly selecting a subset of the data set
        % % Fitting a model to the selected subset
        % % Determining the number of outliers
        % % Repeating steps 1-3 for a prescribed number of iterations
        % % For example, the equation of a line that best fits a set of points can be estimated using RANSAC.
        %
        %     modelLeastSquare = polyfit(tmpData(tmpProstate>0.5),tmpReference(tmpProstate>0.5),1)
        %     points = [tmpData(tmpProstate>0.5),tmpReference(tmpProstate>0.5)];
        %
        %     sampleSize = int64(0.75*length(tmpData(tmpProstate>0.5))); % very random right now
        %     maxDistance = mean(sum((points(:,2)-polyval(modelLeastSquare,points(:,1))).^2,2)) + 2*std(sum((points(:,2)-polyval(modelLeastSquare,points(:,1))).^2,2));
        %
        %     fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
        %     evalLineFcn = @(modelLeastSquare, points) sum((points(:, 2) - polyval(modelLeastSquare, points(:,1))).^2,2);
        %
        %     [modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn, sampleSize,maxDistance);
        %     modelInliers = polyfit(points(inlierIdx,1),points(inlierIdx,2),1)
        %     inlierPts = points(inlierIdx,:);
        %     x = [min(inlierPts(:,1)) max(inlierPts(:,1))];
        %     y = modelInliers(1)*x + modelInliers(2);
        %
        %     figure
        %     plot(inlierPts(:,1),inlierPts(:,2),'*c')
        %     hold on
        %     plot(mdl)
        %     plot(x, y, 'g-')
        %     plot(tmpData(tmpProstate>0.5),tmpData(tmpProstate>0.5)*breg,':m')
        %
        %     title(['c' num2str(cvali) '-map'])
        
        
        R2(cvali,repi)= stats(1);
        f(cvali,repi) = breg(1);
        f_intercept(cvali,repi) = mdl.Coefficients.Estimate(2);
        intercept(cvali, repi) = mdl.Coefficients.Estimate(1);
        dataTilde(:,:,:,cvali,:)= dataM(:,:,:,cvali,:)*breg(1);
        
    end
end

R2out = median(R2,2)
R2iqr = iqr(R2,2)

fout = median(f,2)
fout_intercept = median(f_intercept, 2)
intercept_out = median(intercept,2)
fiqr = iqr(f,2)

end

