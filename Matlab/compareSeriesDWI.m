function [MAD]= compareSeriesDWI(reference, prediction,  scalingFactor, createFigures)
%Karoline Kallis 4/14/2023
% function to compare images series to each other

% Mean absolute difference = MAD
% Mean cross correlation = MX
% MAD = zeros(1,numberOfCompartment);
% MX = zeros(1,numberOfCompartment);


rsi_predicted = zeros([size(prediction.data.outmaps{1}) length(prediction.data.outmaps)]);

rsi = zeros([size(reference.data.outmaps{1}) length(reference.data.outmaps)]);

if(size(rsi_predicted)~=size(rsi))
    disp('Error: Matrices have a different size!!!')
end

ProstateSegmentation = zeros([size(reference.data.outmaps{1}) length(reference.data.outmaps)]);

for pati=1:length(prediction.data.outmaps)
    for cvali=1:size(rsi_predicted,4)
        
        rsi_predicted(:,:,:,cvali,pati) = scalingFactor(cvali)*prediction.data.outmaps_normalized{pati}(:,:,:,cvali);
        rsi(:,:,:,cvali,pati) = reference.data.outmaps_normalized{pati}(:,:,:,cvali);
        %     pati
        ProstateSegmentation(:,:,:,cvali,pati) = reference.data.outseg{pati};
    end
    
    if((sum(rsi_predicted(:,:,:,:,pati),'all')==0 || sum(rsi(:,:,:,:,pati),'all')==0 || sum(ProstateSegmentation(:,:,:,:,pati),'all')==0))
        disp('Dataset empty for patient: ', pati)
    end
    
end



%% Calcualte some values
MAD = mean(abs(rsi_predicted(ProstateSegmentation>0.5)-rsi(ProstateSegmentation>0.5))./rsi(ProstateSegmentation>0.5),'all');

%% Plot images

if(createFigures)
    %     for pati=1:20:size(dwi,5)
    for pati= [4,9,12,23]
        %     for slicei = 10:2:size(CmapsReference,3)-10
        compareDWI(rsi(:,:,:,:,pati).*ProstateSegmentation(:,:,:,:,pati),rsi_predicted(:,:,:,:,pati).*ProstateSegmentation(:,:,:,:,pati), 17, pati)
        %     end
    end
end

%% Plot Boxplots
maxrsi = zeros(size(rsi,5),size(rsi,4));
maxrsi_predicted = zeros(size(rsi,5),size(rsi,4));

minrsi = zeros(size(rsi,5),size(rsi,4));
minrsi_predicted = zeros(size(rsi,5),size(rsi,4));

meanrsi = zeros(size(rsi,5),size(rsi,4));
meanrsi_predicted = zeros(size(rsi,5),size(rsi,4));

medianrsi = zeros(size(rsi,5),size(rsi,4));
medianrsi_predicted = zeros(size(rsi,5),size(rsi,4));

MAD_pati = zeros(size(rsi,5),size(rsi,4));


for bvali = 1:size(rsi,4)
    for pati = 1:size(rsi,5)
        tmprsi = rsi(:,:,:,bvali,pati);
        tmprsiPredicted = rsi_predicted(:,:,:,bvali,pati);
        
        maxrsi(pati,bvali) =  max(tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        maxrsi_predicted(pati,bvali) = max(tmprsiPredicted(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        
        minrsi(pati,bvali) =  min(tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        minrsi_predicted(pati,bvali) = min(tmprsiPredicted(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        
        meanrsi(pati,bvali) =  mean(tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        meanrsi_predicted(pati,bvali) = mean(tmprsiPredicted(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        
        medianrsi(pati,bvali) =  median(tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        medianrsi_predicted(pati,bvali) = median(tmprsiPredicted(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        
        MAD_pati(pati,bvali) = mean(abs(tmprsiPredicted(ProstateSegmentation(:,:,:,bvali,pati)>0.5)-tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5)),'all');
        
    end
end
%figure
%boxplot([maxrsi,maxrsi_predicted, maxrsi-maxrsi_predicted])
if(createFigures)
    
    colors = [0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840];
    grp = [repmat(["TEmin"], 1, length(maxrsi)) repmat(["TEmin_predicted"], 1, length(maxrsi))...
        repmat(["Difference"], 1, length(maxrsi))];
    
    fig = figure;
    %C = [maxC1Prostate{1}; maxC1Prostate{2} ;maxC1Prostate{3} ;maxC1Prostate{4}];
    for bvali = 1:size(rsi,4)
        subplot(1,size(rsi,4),bvali)
        
        boxplot([maxrsi(:,bvali);maxrsi_predicted(:,bvali);maxrsi(:,bvali)-maxrsi_predicted(:,bvali)],grp)
        ylabel('Maximum SI')
        ylim([-1,200])
        title(["c-map: " bvali])
        
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 4]);
    
    %print(fig, 'images/Box_Maximum_TE2min_TE90predicted_cancer','-dpng','-r300');
    
    % fig = figure;
    % %C = [maxC1Prostate{1}; maxC1Prostate{2} ;maxC1Prostate{3} ;maxC1Prostate{4}];
    % for bvali = 1:size(rsi,4)
    %     subplot(1,size(rsi,4),bvali)
    %
    %     boxplot([minrsi(:,bvali);minrsi_predicted(:,bvali);minrsi(:,bvali)-minrsi_predicted(:,bvali)],grp)
    %     ylabel('Minimum SI')
    %     ylim([-2,7])
    %     title(["b-value: " bvali])
    %
    %     h = findobj(gca,'Tag','Box');
    %     for j=1:length(h)
    %         patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    %     end
    % end
    
    fig = figure;
    %C = [maxC1Prostate{1}; maxC1Prostate{2} ;maxC1Prostate{3} ;maxC1Prostate{4}];
    for bvali = 1:size(rsi,4)
        subplot(1,size(rsi,4),bvali)
        
        boxplot([meanrsi(:,bvali);meanrsi_predicted(:,bvali);meanrsi(:,bvali)-meanrsi_predicted(:,bvali)],grp)
        ylabel('Mean SI')
        ylim([-1,30])
        title(["C-map: " bvali])
        
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 4]);
    
    %print(fig, 'images/Box_Mean_TE2min_TE90predicted_cancer','-dpng','-r300');
    % fig = figure;
    % %C = [maxC1Prostate{1}; maxC1Prostate{2} ;maxC1Prostate{3} ;maxC1Prostate{4}];
    % for bvali = 1:size(rsi,4)
    %     subplot(1,size(rsi,4),bvali)
    %
    %     boxplot([medianrsi(:,bvali);medianrsi_predicted(:,bvali);medianrsi(:,bvali)-medianrsi_predicted(:,bvali)],grp)
    %     ylabel('Median SI')
    %     ylim([-2,7])
    %     title(["b-value: " bvali])
    %
    %     h = findobj(gca,'Tag','Box');
    %     for j=1:length(h)
    %         patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    %     end
    % end
    
    
    fig = figure;
    % grpBval = [repmat(["b-value = 0"], 1, length(MAD_pati));...
    %     repmat(["b-value = 50"], 1, length(MAD_pati));...
    %     repmat(["b-value = 800"], 1, length(MAD_pati));...
    %     repmat(["b-value = 1500"], 1, length(MAD_pati));...
    %     repmat(["b-value = 3000"], 1, length(MAD_pati))]';
    
    boxplot(MAD_pati,'Labels',{'c1-map','c2-map','c3-map','c4-map'})
    
    ylabel('MAD [%]')
    %ylim([0,1])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
    
    %print(fig, 'images/Box_MAD_TE2min_TE90predicted_cancer','-dpng','-r300');
    %title(["b-value: " bvali])
end
%% Maximum
mu = (maxrsi+maxrsi_predicted)/2;
varMax =1/2*((maxrsi-mu).^2+ (maxrsi_predicted-mu).^2); % For only 2 oberservations


wSDMax =sqrt( mean(varMax));
RCMax = 2.77*wSDMax

wCVMax= sqrt(mean(varMax./(mu.^2)));
perRCMax = 2.77*wCVMax


%% Compare ADC

% bvals = [0;50;800;1500;3000];
% ADC1 = zeros([size(rsi,1) size(rsi,2) size(rsi,3) size(rsi,5)]);
% ADC2 = zeros([size(rsi,1) size(rsi,2) size(rsi,3) size(rsi,5)]);
% 
% for pati=1:size(rsi,5)
%     ADC1(:,:,:,pati) = compute_ADCs(rsi(:,:,:,:,pati), bvals);
%     ADC2(:,:,:,pati) = compute_ADCs(rsi_predicted(:,:,:,:,pati), bvals);
%     
% end
% 
if(createFigures)
%     
%     %for pati=1:20:size(rsi,5)
%     for pati=[4,9,12,23]
%         
%         %     for slicei = 10:2:size(CmapsReference,3)-10
%         comparersi(ADC1(:,:,:,pati),ADC2(:,:,:,pati), 17, pati)
%     end
%     
%     
    %% Bland altman plot
    for bvali=1:size(maxrsi,2)
        [rpc, fig, stats] = BlandAltman(maxrsi(:,bvali),maxrsi_predicted(:,bvali));
        %print(fig, 'images/BlandAltman_TE2min_TE2min','-dpng','-r300');
    end
end
end

