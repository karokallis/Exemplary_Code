function [MAD]= compareSeriesDWI(prediction, reference, rep_reference,  scalingFactor, createFigures, plotROC)
%Karoline Kallis 4/14/2023
% function to compare images series to each other
% prediction = TE90 series
% Mean absolute difference = MAD
% Mean cross correlation = MX
% MAD = zeros(1,numberOfCompartment);
% MX = zeros(1,numberOfCompartment);

datapath = '/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/PDS_MRI/proc/';

rsi_predicted = zeros([size(reference.data.outmaps{1}) length(reference.data.outmaps)]);
rsi_TE90 = zeros([size(reference.data.outmaps{1}) length(reference.data.outmaps)]);

rsi = zeros([size(reference.data.outmaps{1}) length(reference.data.outmaps)]);
rsi_rep = zeros([size(reference.data.outmaps{1}) length(reference.data.outmaps)]);

if(size(rsi_predicted)~=size(rsi))
    disp('Error: Matrices have a different size!!!')
end

ProstateSegmentation = zeros([size(reference.data.outmaps{1}) length(reference.data.outmaps)]);

for pati=1:length(prediction.data.outmaps)
    for cvali=1:size(rsi_predicted,4)
%         size(prediction.data.outmaps{pati})
%         size(reference.data.outmaps{pati})
        if any(size(prediction.data.outmaps{pati})~=size(reference.data.outmaps{pati}))       
%             continue;
            rsi_predicted(:,:,:,cvali,pati) = imresize3(scalingFactor(cvali)*prediction.data.outmaps_normalized{pati}(:,:,:,cvali), size(reference.data.outmaps_normalized{pati}(:,:,:,cvali)));
            rsi_TE90(:,:,:,cvali,pati) = imresize3(prediction.data.outmaps_normalized{pati}(:,:,:,cvali), size(reference.data.outmaps_normalized{pati}(:,:,:,cvali)));
            %rsi_TE90(:,:,:,cvali,pati) =  vol_resample(prediction.data.outmaps_normalized{pati}(:,:,:,cvali),reference.data.outmaps_normalized{pati}(:,:,:,cvali),eye(4));
        else
            rsi_predicted(:,:,:,cvali,pati) = scalingFactor(cvali)*prediction.data.outmaps_normalized{pati}(:,:,:,cvali);
            rsi_TE90(:,:,:,cvali,pati) =  prediction.data.outmaps_normalized{pati}(:,:,:,cvali);
        end
        %rsi_predicted(:,:,:,cvali,pati) = scalingFactor(cvali)*prediction.data.outmaps_normalized{pati}(:,:,:,cvali);
        rsi(:,:,:,cvali,pati) = reference.data.outmaps_normalized{pati}(:,:,:,cvali);
        rsi_rep(:,:,:,cvali,pati) = rep_reference.data.outmaps_normalized{pati}(:,:,:,cvali);
        ProstateSegmentation(:,:,:,cvali,pati) = reference.data.outseg{pati};
    end
    
    if((sum(rsi_predicted(:,:,:,:,pati),'all')==0 || sum(rsi(:,:,:,:,pati),'all')==0 || sum(ProstateSegmentation(:,:,:,:,pati),'all')==0))
        disp('Dataset empty for patient: ', pati)
    end
    
end


%% Histogram
if createFigures
    colors = [0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840;0.5 0.5 0.5;0.5 0.5 0.0];

    fig1 = figure;
    for cvali=1:size(rsi_predicted,4)
        subplot(2,2,cvali)
        tmp = squeeze(rsi_TE90(:,:,:,cvali,:));
        seg = squeeze(ProstateSegmentation(:,:,:,cvali,:));

        histogram(tmp(seg>0.5 & tmp > 0.0),500,'FaceColor',colors(cvali,:),'EdgeColor',colors(cvali,:))
        if cvali ==1
            xlim([-1,15])
        else
            xlim([-10,150])
        end
        ylim([0, 25000])
    end
    xlabel('SI within prostate')
    xlabel('Frequency')
    sgtitle('TE90')



    fig2 = figure;
    for cvali=1:size(rsi_predicted,4)
        subplot(2,2,cvali)
        tmp = squeeze(rsi_predicted(:,:,:,cvali,:));
        seg = squeeze(ProstateSegmentation(:,:,:,cvali,:));

        histogram(tmp(seg>0.5 & tmp > 0.0),500,'FaceColor',colors(cvali,:),'EdgeColor',colors(cvali,:))
        if cvali ==1
            xlim([-1,15])
        else
            xlim([-10,150])
        end
        ylim([0, 25000])
    end
    xlabel('SI within prostate')
    xlabel('Frequency')
    sgtitle('TEmin*')


    fig3 = figure;
    for cvali=1:size(rsi_predicted,4)
        subplot(2,2,cvali)
        tmp = squeeze(rsi(:,:,:,cvali,:));
        seg = squeeze(ProstateSegmentation(:,:,:,cvali,:));

        histogram(tmp(seg>0.5 & tmp > 0.0),500,'FaceColor',colors(cvali,:),'EdgeColor',colors(cvali,:))
        if cvali ==1
            xlim([-1,15])
        else
            xlim([-10,150])
        end
        ylim([0, 25000])
    end
    xlabel('SI within prostate')
    xlabel('Frequency')
    sgtitle('TEmin1')


    fig4 = figure;
    for cvali=1:size(rsi_predicted,4)
        subplot(2,2,cvali)
        tmp = squeeze(rsi_rep(:,:,:,cvali,:));
        seg = squeeze(ProstateSegmentation(:,:,:,cvali,:));

        histogram(tmp(seg>0.5 & tmp > 0.0),500,'FaceColor',colors(cvali,:),'EdgeColor',colors(cvali,:))
        if cvali ==1
            xlim([-1,15])
        else
            xlim([-10,150])
        end
        ylim([0, 25000])
    end
    xlabel('SI within prostate')
    xlabel('Frequency')
    sgtitle('TEmin2')
end
%% Calcualte some values
% 1 = TE90_cali - RSI
% 2 = TE90_cali - RSI_REP
% 3 = RSI_REP - RSI
% 4 = TE90 - RSI
% 5 = TE90 - TE90_calibrated
MAD = zeros(1,5);
MAD(1) = mean(abs(rsi_predicted(ProstateSegmentation>0.5)-rsi(ProstateSegmentation>0.5)),'all');
MAD(2) = mean(abs(rsi_predicted(ProstateSegmentation>0.5)-rsi_rep(ProstateSegmentation>0.5)),'all');
MAD(3) = mean(abs(rsi_rep(ProstateSegmentation>0.5)-rsi(ProstateSegmentation>0.5)),'all');
MAD(4) = mean(abs(rsi_TE90(ProstateSegmentation>0.5)-rsi(ProstateSegmentation>0.5)),'all');
MAD(5) = mean(abs(rsi_TE90(ProstateSegmentation>0.5)-rsi_predicted(ProstateSegmentation>0.5)),'all');

%% Plot images

% if(createFigures)
%     %for pati=1:1:size(rsi,5)
%     for pati= [4,9,12,23]
%         %     for slicei = 10:2:size(CmapsReference,3)-10
% %         pati
% %         compareDWI(rsi(:,:,:,:,pati),rsi_TE90(:,:,:,:,pati), 17, pati,'TE')
% %         pause
%         compareDWI(rsi(:,:,:,:,pati).*ProstateSegmentation(:,:,:,:,pati),rsi_predicted(:,:,:,:,pati).*ProstateSegmentation(:,:,:,:,pati), 17, pati,'TEminvsTE90*')
%         compareDWI(rsi(:,:,:,:,pati).*ProstateSegmentation(:,:,:,:,pati),rsi_TE90(:,:,:,:,pati).*ProstateSegmentation(:,:,:,:,pati), 17, pati,'TEminvsTE90')
%         compareDWI(rsi_predicted(:,:,:,:,pati).*ProstateSegmentation(:,:,:,:,pati),rsi_TE90(:,:,:,:,pati).*ProstateSegmentation(:,:,:,:,pati), 17, pati,'TE90*vsTE90')
%
%         %     end
%     end
% end

%% Plot Boxplots
maxrsi = zeros(size(rsi,5),size(rsi,4));
maxrsi_predicted = zeros(size(rsi,5),size(rsi,4));
maxrsi_rep = zeros(size(rsi,5),size(rsi,4));
maxrsi_TE90 = zeros(size(rsi,5),size(rsi,4));

meanrsi = zeros(size(rsi,5),size(rsi,4));
meanrsi_predicted = zeros(size(rsi,5),size(rsi,4));
meanrsi_rep = zeros(size(rsi,5),size(rsi,4));
meanrsi_TE90 = zeros(size(rsi,5),size(rsi,4));

rsi_98 = zeros(size(rsi,5),size(rsi,4));
rsi_98_predicted = zeros(size(rsi,5),size(rsi,4));
rsi_98_rep = zeros(size(rsi,5),size(rsi,4));
rsi_98_TE90 = zeros(size(rsi,5),size(rsi,4));

% MAD_pati_rsi_pred = zeros(size(rsi,5),size(rsi,4));
% MAD_pati_rsi_rep = zeros(size(rsi,5),size(rsi,4));
% MAD_pati_rsi_TE90 = zeros(size(rsi,5),size(rsi,4));
% MAD_pati_rep_TE90 = zeros(size(rsi,5),size(rsi,4));
% MAD_pati_rep_pred = zeros(size(rsi,5),size(rsi,4));



for bvali = 1:size(rsi,4)
    for pati = 1:size(rsi,5)
        tmprsi = rsi(:,:,:,bvali,pati);
        tmprsiPredicted = rsi_predicted(:,:,:,bvali,pati);
        tmprsi_rep = rsi_rep(:,:,:,bvali,pati);
        tmprsi_TE90 = rsi_TE90(:,:,:,bvali,pati);
        
        
        maxrsi(pati,bvali) =  max(tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        maxrsi_predicted(pati,bvali) = max(tmprsiPredicted(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        maxrsi_rep(pati,bvali) =  max(tmprsi_rep(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        maxrsi_TE90(pati,bvali) = max(tmprsi_TE90(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        
        meanrsi(pati,bvali) =  mean(tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        meanrsi_predicted(pati,bvali) = mean(tmprsiPredicted(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        meanrsi_rep(pati,bvali) =  mean(tmprsi_rep(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        meanrsi_TE90(pati,bvali) = mean(tmprsi_TE90(ProstateSegmentation(:,:,:,bvali,pati)>0.5));
        
        rsi_98(pati,bvali) =  prctile(tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5),98);
        rsi_98_predicted(pati,bvali) = prctile(tmprsiPredicted(ProstateSegmentation(:,:,:,bvali,pati)>0.5),98);
        rsi_98_rep(pati,bvali) =  prctile(tmprsi_rep(ProstateSegmentation(:,:,:,bvali,pati)>0.5),98);
        rsi_98_TE90(pati,bvali) = prctile(tmprsi_TE90(ProstateSegmentation(:,:,:,bvali,pati)>0.5),98);
        %
        %         MAD_pati_rsi_pred(pati,bvali) = mean(abs(tmprsiPredicted(ProstateSegmentation(:,:,:,bvali,pati)>0.5)-tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5)),'all');
        %         MAD_pati_rsi_rep(pati,bvali) = mean(abs(tmprsi_rep(ProstateSegmentation(:,:,:,bvali,pati)>0.5)-tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5)),'all');
        %         MAD_pati_rsi_TE90(pati,bvali) = mean(abs(tmprsi_TE90(ProstateSegmentation(:,:,:,bvali,pati)>0.5)-tmprsi(ProstateSegmentation(:,:,:,bvali,pati)>0.5)),'all');
        %         MAD_pati_rep_TE90(pati,bvali) = mean(abs(tmprsi_TE90(ProstateSegmentation(:,:,:,bvali,pati)>0.5)-tmprsi_rep(ProstateSegmentation(:,:,:,bvali,pati)>0.5)),'all');
        %         MAD_pati_rep_pred(pati,bvali) = mean(abs(tmprsiPredicted(ProstateSegmentation(:,:,:,bvali,pati)>0.5)-tmprsi_rep(ProstateSegmentation(:,:,:,bvali,pati)>0.5)),'all');
        
    end
end
%figure
%boxplot([maxrsi,maxrsi_predicted, maxrsi-maxrsi_predicted])
if(createFigures)
    
    colors = [0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840;0.5 0.5 0.5;0.5 0.5 0.0];
    grp = [repmat(["TE90"], 1, length(maxrsi)) repmat(["TEmin*"], 1, length(maxrsi))...
        repmat(["TEmin"], 1, length(maxrsi)), repmat(["TEmin_rep"], 1, length(maxrsi))];
%     
    fig = figure;
    %C = [maxC1Prostate{1}; maxC1Prostate{2} ;maxC1Prostate{3} ;maxC1Prostate{4}];
    for bvali = 1:size(rsi,4)
        subplot(1,size(rsi,4),bvali)
        %violin([maxrsi_TE90(:,bvali),maxrsi_predicted(:,bvali), maxrsi(:,bvali),maxrsi_rep(:,bvali)],[],'facecolor',repmat(colors,5,1),'mc',[],'medc','k');
        
        boxplot([maxrsi_TE90(:,bvali);maxrsi_predicted(:,bvali); maxrsi(:,bvali);maxrsi_rep(:,bvali)],grp)
        ylabel('Maximum SI')
        if bvali==1
            ylim([-1,100])
        else
            ylim([-1,300])
        end
        title(["c-map: " bvali])
        
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 4]);  
    print(fig, 'images/Box_Maximum_RSI_Focus','-dpng','-r300');
    
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
        %violin([meanrsi_TE90(:,bvali),meanrsi_predicted(:,bvali), meanrsi(:,bvali),meanrsi_rep(:,bvali)],[],'facecolor',repmat(colors,5,1),'mc',[],'medc','k');
        
        boxplot([meanrsi_TE90(:,bvali);meanrsi_predicted(:,bvali); meanrsi(:,bvali);meanrsi_rep(:,bvali)],grp)
        ylabel('Mean SI')
        if bvali==1
            ylim([-1,20])
        else
            ylim([-1,100])
        end
        title(["C-map: " bvali])
        
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 4]);
    print(fig, 'images/Box_Mean_RSI_Focus','-dpng','-r300');
    
    
    fig = figure;
    %C = [maxC1Prostate{1}; maxC1Prostate{2} ;maxC1Prostate{3} ;maxC1Prostate{4}];
    for bvali = 1:size(rsi,4)
        subplot(1,size(rsi,4),bvali)
        %violin([rsi_98_TE90(:,bvali),rsi_98_predicted(:,bvali), rsi_98(:,bvali),rsi_98_rep(:,bvali)],[],'facecolor',repmat(colors,5,1),'mc',[],'medc','k');
        
        boxplot([rsi_98_TE90(:,bvali);rsi_98_predicted(:,bvali); rsi_98(:,bvali);rsi_98_rep(:,bvali)],grp)
        ylabel('98th percentile SI')
        if bvali==1
            ylim([-1,100])
        else
            ylim([-1,300])
        end
        title(["C-map: " bvali])
        
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 4]);
    print(fig, 'images/Box_98thPer_RSI_Focus_cancer','-dpng','-r300');
    
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
    
    %
    %     fig = figure;
    %     %C = [maxC1Prostate{1}; maxC1Prostate{2} ;maxC1Prostate{3} ;maxC1Prostate{4}];
    %     for bvali = 1:size(rsi,4)
    %         subplot(1,size(rsi,4),bvali)
    %     % grpBval = [repmat(["b-value = 0"], 1, length(MAD_pati));...
    %     %     repmat(["b-value = 50"], 1, length(MAD_pati));...
    %     %     repmat(["b-value = 800"], 1, length(MAD_pati));...
    %     %     repmat(["b-value = 1500"], 1, length(MAD_pati));...
    %     %     repmat(["b-value = 3000"], 1, length(MAD_pati))]';
    %
    %     boxplot(MAD_pati_rsi_rep,'Labels',{'c1-map','c2-map','c3-map','c4-map'})
    %
    %     ylabel('MAD SI')
    %     %ylim([0,1])
    %     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
    
    %print(fig, 'images/Box_MAD_TE2min_TE90predicted_cancer','-dpng','-r300');
    %title(["b-value: " bvali])
end

%% Histograms
if(createFigures)
    
    fig = figure;
    %C = [maxC1Prostate{1}; maxC1Prostate{2} ;maxC1Prostate{3} ;maxC1Prostate{4}];
    
    
    for bvali = 1:size(rsi,4)
        if bvali==1
            x = 0:0.01:20;
        else
            x = 0:0.01:300;
        end
        subplot(1,size(rsi,4),bvali)
        pd = fitdist(rsi_98_TE90(:,bvali),'kernel');
        y =pdf(pd,x);
        plot(x, y, 'Linewidth',2)
        hold on
        pd = fitdist(rsi_98_predicted(:,bvali),'kernel');
        y =pdf(pd,x);
        plot(x, y, 'Linewidth',2)
        
        pd = fitdist(rsi_98(:,bvali),'kernel');
        y =pdf(pd,x);
        plot(x, y, 'Linewidth',2)
        
        pd = fitdist(rsi_98_rep(:,bvali),'kernel');
        y =pdf(pd,x);
        plot(x, y, 'Linewidth',2)
        title(["C-map: " bvali])
        ylabel('PDF')
        xlabel('98th within prostate')
    end
    legend('TE90','TEmin*','TEmin_1', 'TEmin_2')
    
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 4]);  
    print(fig, 'images/Hist_98th_RSI_Focus_cancer','-dpng','-r300');

    
    fig = figure;
    for bvali = 1:size(rsi,4)
        if bvali==1
            x = 0:0.01:5;
        else
            x = 0:0.01:100;
        end
        subplot(1,size(rsi,4),bvali)
        pd = fitdist(meanrsi_TE90(:,bvali),'kernel');
        y =pdf(pd,x)/100*(max(meanrsi_TE90(:,bvali)) - min(meanrsi_TE90(:,bvali)));
        plot(x, y, 'Linewidth',2)
        hold on
        pd = fitdist(meanrsi_predicted(:,bvali),'kernel');
        y =pdf(pd,x)/100*(max(meanrsi_predicted(:,bvali)) - min(meanrsi_predicted(:,bvali)));
        plot(x, y, 'Linewidth',2)
        
        pd = fitdist(meanrsi(:,bvali),'kernel');
        y =pdf(pd,x)/100*(max(meanrsi(:,bvali)) - min(meanrsi(:,bvali)));
        plot(x, y, 'Linewidth',2)
        
        pd = fitdist(meanrsi_rep(:,bvali),'kernel');
        y =pdf(pd,x)/100*(max(meanrsi_rep(:,bvali)) - min(meanrsi_rep(:,bvali)));
        plot(x, y, 'Linewidth',2)
        title(["C-map: " bvali])
        ylabel('PDF')
        xlabel('Mean within prostate')
    end
    legend('TE90','TEmin*','TEmin_1', 'TEmin_2')
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 4]);  
    print(fig, 'images/Hist_Mean_RSI_Focus_cancer','-dpng','-r300');
    
    fig = figure;
    for bvali = 1:size(rsi,4)
        if bvali==1
            x = 0:0.01:50;
        else
            x = 0:0.01:400;
        end
        subplot(1,size(rsi,4),bvali)
        pd = fitdist(maxrsi_TE90(:,bvali),'kernel');
        y =pdf(pd,x);
        plot(x, y, 'Linewidth',2)
        hold on
        pd = fitdist(maxrsi_predicted(:,bvali),'kernel');
        y =pdf(pd,x);
        plot(x, y, 'Linewidth',2)
        
        pd = fitdist(maxrsi(:,bvali),'kernel');
        y =pdf(pd,x);
        plot(x, y, 'Linewidth',2)
        
        pd = fitdist(maxrsi_rep(:,bvali),'kernel');
        y =pdf(pd,x);
        plot(x, y, 'Linewidth',2)
        title(["C-map: " bvali])
        ylabel('PDF')
        xlabel('Maximum within prostate')
    end
    legend('TE90','TEmin*','TEmin_1', 'TEmin_2')
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 4]);  
    print(fig, 'images/Hist_Maximum_RSI_Focus_cancer','-dpng','-r300');
end

%% ROC

figName = 'ROC Maximum';
fignumber = 6666;
% downsample data


if plotROC
    csPCa = double(cell2mat(reference.data.csPCa_status));
    idx = [randsample(find(csPCa==0),sum(csPCa)) find(csPCa==1)];
%     Ind = f(randi(length(f)));
%     benignIdx = find(csPCa 
    
    predictor{1} = maxrsi_TE90(idx,1);
    predictor{2} = maxrsi_predicted(idx,1);
    predictor{3} = maxrsi(idx,1);
    predictor{4} = maxrsi_rep(idx,1);

%     predictor{1} = rsi_98_TE90(idx,1);
%     predictor{2} = rsi_98_predicted(idx,1);
%     predictor{3} = rsi_98(idx,1);
%     predictor{4} = rsi_98_rep(idx,1);
    
    
    % predictor{length(datafiles)+1} = cellfun(@(x,y) max(x.*y,[],'all','omitnan'), data.ADCmap, data.outseg, 'UniformOutput', true)';
    
    for pp=1:length(predictor)
        
        tmproc = tms_plot_ROC(csPCa(idx),predictor{pp},'rocplotflag',0,'XVals','all','UseNearest','on','reverseflag',0,'reversewarning',0);
        
        xx{pp}=tmproc.x;
        yy{pp}=tmproc.y;
        auc{pp}=tmproc.AUC;
        FPR90{pp}=tmproc.FPR90;
    end
    % legendnames = {'RSI -- AUC=0.70','aDWI -- AUC=0.59','sDWI_{1000} -- AUC=0.45','sDWI_{500} -- AUC=0.47'};
    legendnames = {'TE90','TEmin*','TEmin_1','TEmin_2'};
    
    tms_plot_ROC_curves(xx,yy,'outfile',['images/ROC',figName],'saveflag',1,'legendnames',legendnames,...
        'titlestr',figName,'fignum',fignumber,'forceflag',0 ,...
        'linestylelist',{'-','-','-','-'},'markerlist',{'none','none','none','none'},...
        'colorlist',{[1 0.0 0.0],[0 0 0],[0 0 1],[0.0 0.4 0.4]});
    
end

% %% Maximum
% mu = (maxrsi+maxrsi_predicted)/2;
% varMax =1/2*((maxrsi-mu).^2+ (maxrsi_predicted-mu).^2); % For only 2 oberservations
% 
% 
% wSDMax =sqrt( mean(varMax));
% RCMax = 2.77*wSDMax
% 
% wCVMax= sqrt(mean(varMax./(mu.^2)));
% perRCMax = 2.77*wCVMax

%% Combo Plot -- Maximum 

fig5 = figure;

colors = [0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840;0.5 0.5 0.5;0.5 0.5 0.0];
grp = [repmat(["TE90"], 1, length(maxrsi)) repmat(["TEmin*"], 1, length(maxrsi))...
    repmat(["TEmin1"], 1, length(maxrsi)), repmat(["TEmin2"], 1, length(maxrsi))];

for bvali = 1:size(rsi,4)
    subplot(size(rsi,4),2,2*bvali-1)
    boxplot([maxrsi_TE90(:,bvali);maxrsi_predicted(:,bvali); maxrsi(:,bvali);maxrsi_rep(:,bvali)],grp)
    ylabel('Maximum SI')
    if bvali==1
        ylim([-1,100])
    else
        ylim([-1,300])
    end
    %title(["c-map: " bvali])

    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.8);
    end
    %     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 4]);
    %     print(fig, 'images/Box_Maximum_RSI_Focus_cancer','-dpng','-r300');

    subplot(size(rsi,4),2,2*bvali)
    if bvali==1
        x = 0:0.01:100;
    else
        x = 0:0.01:300;
        ylim([-1,300])
    end

    pd = fitdist(maxrsi_TE90(:,bvali),'kernel');
    y =pdf(pd,x);
    %y =pdf(pd,x)/100*(max(maxrsi_TE90(:,bvali)) - min(maxrsi_TE90(:,bvali)));
    plot(y, x, 'Linewidth',2, 'Color',[0.5 0.5 0.5])
    hold on
    pd = fitdist(maxrsi_predicted(:,bvali),'kernel');
    y =pdf(pd,x);
    %y =pdf(pd,x)/100*(max(maxrsi_predicted(:,bvali)) - min(maxrsi_predicted(:,bvali)));
    plot(y, x, 'Linewidth',2,'Color',[0.6350 0.0780 0.1840])

    pd = fitdist(maxrsi(:,bvali),'kernel');
    y =pdf(pd,x);
    %y =pdf(pd,x)/100*(max(maxrsi(:,bvali)) - min(maxrsi(:,bvali)));
    plot(y, x, 'Linewidth',2,'Color',[0 0.4470 0.7410])

    pd = fitdist(maxrsi_rep(:,bvali),'kernel');
    y =pdf(pd,x);
    %y =pdf(pd,x)/100*(max(maxrsi_rep(:,bvali)) - min(maxrsi_rep(:,bvali)));
    plot(y, x, 'Linewidth',2,'Color',[0.9290 0.6940 0.1250])
    %title(["c-map: " bvali])
    %ylabel('PDF')
    %xlabel('Maximum SI')
    %xlim([0,0.08])
    if bvali==1
        ylim([-1,100])
        legend('TE90','TEmin*','TEmin1', 'TEmin2')
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 9]);  
    print(fig5, 'images/Fig98th_RSI_Focus','-dpng','-r600');

end


%% Combo Plot -- 98th Percentile 

fig6 = figure;

colors = [0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840;0.5 0.5 0.5;0.5 0.5 0.0];
grp = [repmat(["TE90"], 1, length(maxrsi)) repmat(["TEmin*"], 1, length(maxrsi))...
    repmat(["TEmin1"], 1, length(maxrsi)), repmat(["TEmin2"], 1, length(maxrsi))];

for bvali = 1:size(rsi,4)
    subplot(size(rsi,4),2,2*bvali-1)
    boxplot([rsi_98_TE90(:,bvali);rsi_98_predicted(:,bvali); rsi_98(:,bvali);rsi_98_rep(:,bvali)],grp)
    ylabel('98th percentile SI')
    if bvali==1
        ylim([-1,50])
    else
        ylim([-1,300])
    end
    %title(["c-map: " bvali])

    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.8);
    end
    %     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 4]);
    %     print(fig, 'images/Box_Maximum_RSI_Focus_cancer','-dpng','-r300');

    subplot(size(rsi,4),2,2*bvali)
    if bvali==1
        x = 0:0.01:50;
        ylim([-1,50])
    else
        x = 0:0.01:300;
        ylim([-1,300])
    end

    pd = fitdist(rsi_98_TE90(:,bvali),'kernel');
    y =pdf(pd,x);
    %y =pdf(pd,x)/100*(max(rsi_98_TE90(:,bvali)) - min(rsi_98_TE90(:,bvali)));
    plot(y, x, 'Linewidth',2, 'Color',[0.5 0.5 0.5])
    hold on
    pd = fitdist(rsi_98_predicted(:,bvali),'kernel');
    y =pdf(pd,x);
    %y =pdf(pd,x)/100*(max(rsi_98_predicted(:,bvali)) - min(rsi_98_predicted(:,bvali)));
    plot(y, x, 'Linewidth',2,'Color',[0.6350 0.0780 0.1840])

    pd = fitdist(rsi_98(:,bvali),'kernel');
    y =pdf(pd,x);
    %y =pdf(pd,x)/100*(max(rsi_98(:,bvali)) - min(rsi_98(:,bvali)));
    plot(y, x, 'Linewidth',2,'Color',[0 0.4470 0.7410])

    pd = fitdist(rsi_98_rep(:,bvali),'kernel');
    y =pdf(pd,x);
    %y =pdf(pd,x)/100*(max(rsi_98_rep(:,bvali)) - min(rsi_98_rep(:,bvali)));
    plot(y, x, 'Linewidth',2,'Color',[0.9290 0.6940 0.1250])
    %title(["c-map: " bvali])
    %ylabel('PDF')
    %xlabel('Maximum SI')
    %xlim([0,0.08])
    if bvali==1
        legend('TE90','TEmin*','TEmin1', 'TEmin2')
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 9]);  
    print(fig6, 'images/Fig98th_RSI_Focus','-dpng','-r600');

end
end

