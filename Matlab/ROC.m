function [auc] =  ROC(TE90, TEmin, TEmin_rep, f, fignumber)
% Created by Karoline Kallis 08/09/2022
% Function to compare ROC curves of diffusion maps
addpath(genpath('../../synBvalues/matlab/ROC/'))

% masked == true if prostate mask is considered and not full FOV
% meta = load(metafile); % loads ptlevel metadata
%data=load(datafiles(3).fname); % Load data each time in order to create different margins

%compareDWISlices(data.map{25}.*data.outseg{25},data.syn1000{25}.*data.outseg{25},data.syn500{25}.*data.outseg{25},8, 'Bvalues')
% if ~isempty(margin)
%     data = dilateContour(data, margin); % margin in mm
% end

%% Create predictors

predictor{1} = cellfun(@(x,y) max(x(y>0.5),[],'all','omitnan'), TE90.data.outmaps_normalized(:,:,:,1), TE90.data.outseg, 'UniformOutput', true)';
predictor{2} = cellfun(@(x,y) max(f(1)*x(y>0.5),[],'all','omitnan'),TE90.data.outmaps_normalized(:,:,:,1), TE90.data.outseg, 'UniformOutput', true)';
predictor{3} = cellfun(@(x,y) max(x(y>0.5),[],'all','omitnan'), TEmin.data.outmaps_normalized(:,:,:,1), TE90.data.outseg, 'UniformOutput', true)';
predictor{4} = cellfun(@(x,y) max(x(y>0.5),[],'all','omitnan'), TEmin_rep.data.outmaps_normalized(:,:,:,1), TE90.data.outseg, 'UniformOutput', true)';



% predictor{length(datafiles)+1} = cellfun(@(x,y) max(x.*y,[],'all','omitnan'), data.ADCmap, data.outseg, 'UniformOutput', true)';

for pp=1:length(predictor)
    csPCa = double(cell2mat(datafiles(pp).data.csPCa_status));
    tmproc = tms_plot_ROC(csPCa,predictor{pp},'rocplotflag',0,'XVals','all','UseNearest','on','reverseflag',0,'reversewarning',0);
    
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



% %% Bootstrapping
% %% Compare patient-level ROC curves with bootstrapping
% %close all
% % cand_model_list = {'RSIrs', 'aDWI','sDWI1000','sDWI500'};
% bootstrapN = 1e4;
% subsetidx = meta.subsetidx_all;
% roititlename = 'Prostate';
% subjects = meta.metasubjectIDs;
% nvox =1;
% fignuminc = 1000;
% outdir ='../data/';
% roiname = 'prostate';
% 
% % Desin predictor
% pred.prostate.RSI_FOCUS = predictor{1};
% pred.prostate.RSI_TE90 = predictor{2};
% pred.prostate.RSI_TEmin1 = predictor{3};
% pred.prostate.RSI_TEmin2 = predictor{4};
% 
% 
% 
% % % pred.(roiname).(modelname)(si,1) = median(roivoxels(1:nvox));
% 
% if bootstrapN
%     bootmodels = {'RSI_FOCUS', 'RSI_TE90','RSI_TEmin1','RSI_TEmin2'}; % {'pirads','quadexp','ADC','geADC'}; % no bootstrapping for DWI
%     
%     for mi=1:numel(bootmodels)
%         if mi==1; models_str = bootmodels{mi};
%         else models_str = sprintf('%s_%s',models_str,bootmodels{mi});
%         end
%     end
%     stratifygroups{1} = subjects(meta.csPCa_status==1 & subsetidx); 
%     stratifygroups{2} = subjects(meta.csPCa_status==0 & subsetidx);
%     
%     fname_boot = sprintf('%s/patientlevel_RSIapp_worstnvox%d_%s_roi%s_noreverseAUC_bootstrapN%d.mat',outdir,nvox,models_str,roititlename,bootstrapN);
%     
% %     if plusPIRADSflag
% %         fname_boot = regexprep(fname_boot,'_roi','_plusPIRADS_roi');
% %         for mi=1:numel(bootmodels)
% %             if ~strcmp(bootmodels{mi},'pirads')
% %                 bootmodels{end+1} = sprintf('%s_plus_PIRADS',bootmodels{mi});
% %             end
% %         end
% %     end
%     
%     fignum = fignuminc;
%     for bi=1:numel(bootmodels)
%         if regexp(bootmodels{bi},'ADC')
%             reverseflags(bi)=2; 
%         else
%             reverseflags(bi)=0;
%         end
%     end
%     [g.bootstats,g.bootsubs,g.b]=patientbootstrap(meta.csPCa_status,pred.prostate,subjects,stratifygroups,bootmodels,...
%         'bootstrapN',bootstrapN,'fname_boot',fname_boot,'fignum',fignum,'plotflag',0,'metriclist',{'AUC','FPR90'},...
%         'subsetidx',subsetidx,'titlestr',roititlename,'reverseflags',reverseflags);
%     
% end

%% Results Prostae only
% Summary statistics for N=10000 bootstrap samples:
% 	 RSIrs: bootstrap median AUC=0.7841 [0.7063,0.8555]. 	mean 0.7833 
% 	 aDWI: bootstrap median AUC=0.6222 [0.5307,0.7095]. 	mean 0.6214 
% 	 sDWI1000: bootstrap median AUC=0.6487 [0.5585,0.7335]. 	mean 0.6479 
% 	 sDWI500: bootstrap median AUC=0.6320 [0.5395,0.7190]. 	mean 0.6316 
% 
% bootstrap difference aDWI vs RSIrs for AUC: two-sided p=0
% 
% bootstrap difference sDWI1000 vs RSIrs for AUC: two-sided p=0.0012
% 
% bootstrap difference sDWI500 vs RSIrs for AUC: two-sided p=0.0018
% 
% bootstrap difference sDWI1000 vs aDWI for AUC: two-sided p=0.621
% 
% bootstrap difference sDWI500 vs aDWI for AUC: two-sided p=0.861
% 
% bootstrap difference sDWI500 vs sDWI1000 for AUC: two-sided p=0.544


%% Prostate + 5mm
% Summary statistics for N=10000 bootstrap samples:
% 	 RSIrs: bootstrap median AUC=0.7671 [0.6872,0.8399]. 	mean 0.7660 
% 	 aDWI: bootstrap median AUC=0.6061 [0.5166,0.6928]. 	mean 0.6058 
% 	 sDWI1000: bootstrap median AUC=0.5993 [0.5098,0.6866]. 	mean 0.5995 
% 	 sDWI500: bootstrap median AUC=0.5614 [0.4691,0.6510]. 	mean 0.5608 
% 
% bootstrap difference aDWI vs RSIrs for AUC: two-sided p=0
% 
% bootstrap difference sDWI1000 vs RSIrs for AUC: two-sided p=0.0014
% 
% bootstrap difference sDWI500 vs RSIrs for AUC: two-sided p=0.0002
% 
% bootstrap difference sDWI1000 vs aDWI for AUC: two-sided p=0.922
% 
% bootstrap difference sDWI500 vs aDWI for AUC: two-sided p=0.521
% 
% bootstrap difference sDWI500 vs sDWI1000 for AUC: two-sided p=0.209

%% wFOV
% Summary statistics for N=10000 bootstrap samples:
% 	 RSIrs: bootstrap median AUC=0.7000 [0.6138,0.7821]. 	mean 0.6997 
% 	 aDWI: bootstrap median AUC=0.5889 [0.4970,0.6808]. 	mean 0.5890 
% 	 sDWI1000: bootstrap median AUC=0.4481 [0.3554,0.5417]. 	mean 0.4477 
% 	 sDWI500: bootstrap median AUC=0.4691 [0.3761,0.5614]. 	mean 0.4688 
% 
% bootstrap difference aDWI vs RSIrs for AUC: two-sided p=0.0036
% 
% bootstrap difference sDWI1000 vs RSIrs for AUC: two-sided p=0.0002
% 
% bootstrap difference sDWI500 vs RSIrs for AUC: two-sided p=0.0004
% 
% bootstrap difference sDWI1000 vs aDWI for AUC: two-sided p=0.0206
% 
% bootstrap difference sDWI500 vs aDWI for AUC: two-sided p=0.0456
% 
% bootstrap difference sDWI500 vs sDWI1000 for AUC: two-sided p=0.234

%% 
% compareHistogram(predictor{2}',predictor{3}', 'Diffusion', 'synB1000', 'FPRSI', '', [])
% compareHistogram(predictor{2}',predictor{4}', 'Diffusion', 'synB500', 'FPRSI', '',[])
% %% Plotting 
% fig= figure;
% %z = [predictor{2}; predictor{3}; predictor{4}];
% %g = [zeros(length(statsCancerFPRSI.Min),1); ones(length(statsCancerKOP.Min),1);...
% %    2*ones(length(statsCancerNGP_TEmin11.Min),1);3*ones(length(statsCancerNGP_TEmin7.Min),1);4*ones(length(statsCancerNGP_TE90.Min),1);5*ones(length(statsCancerNGP_Vendor.Min),1);...
% %    6*ones(length(statsCancerDream_TE80.Min),1);7*ones(length(statsCancerDream_TE100.Min),1);8*ones(length(statsCancerDream_Diff.Min),1)];
% boxplot([predictor{2},predictor{3},predictor{4}],'Labels',{'acquired B2000', 'syn B2000 (B1000)', 'syn B2000 (B500)'},'Colors', 'rcb')
% ylim([0, 3000])
% title('Maximum')
% print('images/Boxplot','-dpng','-r300');

%%

% figure
% subplot(2,1,1)
% qqplot(predictor{2},predictor{3})
% title('acquired B2000 vs. syn B2000 (B1000)')
% 
% subplot(2,1,2)
% qqplot(predictor{2},predictor{4})
% title('acquired B2000  vs. syn B2000 (B500)')
% % 
% % subplot(3,1,3)
% % qqplot(predictor{3},predictor{4})
% % title('syn 1000 vs syn 500')
% print('images/qqplot','-dpng','-r300');

% %% Test Significance
% if masked
%     ben1 = cellfun(@(x,y,z) sqrt(mean((nonzeros(x.*y)-nonzeros(z.*y)).^2,'all')), data.map(meta.csPCa_status==0), data.outseg(meta.csPCa_status==0),data.syn1000(meta.csPCa_status==0),'UniformOutput', true);
%     can1 = cellfun(@(x,y,z) sqrt(mean((nonzeros(x.*y)-nonzeros(z.*y)).^2,'all')), data.map(meta.csPCa_status==1), data.outseg(meta.csPCa_status==1),data.syn1000(meta.csPCa_status==1),'UniformOutput', true);
%     
%     [~, p1] = ttest2(ben1, can1);
%     
%     ben2 = cellfun(@(x,y,z) sqrt(mean((nonzeros(x.*y)-nonzeros(z.*y)).^2,'all')), data.map(meta.csPCa_status==0), data.outseg(meta.csPCa_status==0),data.syn500(meta.csPCa_status==0),'UniformOutput', true);
%     can2 = cellfun(@(x,y,z) sqrt(mean((nonzeros(x.*y)-nonzeros(z.*y)).^2,'all')), data.map(meta.csPCa_status==1), data.outseg(meta.csPCa_status==1),data.syn500(meta.csPCa_status==1),'UniformOutput', true);
%     
%     [~, p2] = ttest2(ben2, can2);
%     
%     ben3 = cellfun(@(x,y,z) sqrt(mean((nonzeros(x.*y)-nonzeros(z.*y)).^2,'all')), data.syn1000(meta.csPCa_status==0), data.outseg(meta.csPCa_status==0),data.syn500(meta.csPCa_status==0),'UniformOutput', true);
%     can3 = cellfun(@(x,y,z) sqrt(mean((nonzeros(x.*y)-nonzeros(z.*y)).^2,'all')), data.syn1000(meta.csPCa_status==1), data.outseg(meta.csPCa_status==1),data.syn500(meta.csPCa_status==1),'UniformOutput', true);
%     
%     [~, p3] = ttest2(ben3, can3);
%     
% else
%     ben1 = cellfun(@(x,z) sqrt(mean((nonzeros(x(x~=0))-nonzeros(z(x~=0))).^2,'all')), data.map(meta.csPCa_status==0),data.syn1000(meta.csPCa_status==0),'UniformOutput', true);
%     can1 = cellfun(@(x,z) sqrt(mean((nonzeros(x(x~=0))-nonzeros(z(x~=0))).^2,'all')), data.map(meta.csPCa_status==1),data.syn1000(meta.csPCa_status==1),'UniformOutput', true);
%     
%     [~, p1] = ttest2(ben1, can1);
%     
%     ben2 = cellfun(@(x,z) sqrt(mean((nonzeros(x(x~=0))-nonzeros(z(x~=0))).^2,'all')), data.map(meta.csPCa_status==0),data.syn500(meta.csPCa_status==0),'UniformOutput', true);
%     can2 = cellfun(@(x,z) sqrt(mean((nonzeros(x(x~=0))-nonzeros(z(x~=0))).^2,'all')), data.map(meta.csPCa_status==1),data.syn500(meta.csPCa_status==1),'UniformOutput', true);
%     
%     [~, p2] = ttest2(ben2, can2);
%     
%     ben3= cellfun(@(x,z) sqrt(mean((nonzeros(x(x~=0))-nonzeros(z(x~=0))).^2,'all')), data.syn1000(meta.csPCa_status==0),data.syn500(meta.csPCa_status==0),'UniformOutput', true);
%     can3 = cellfun(@(x,z) sqrt(mean((nonzeros(x(x~=0))-nonzeros(z(x~=0))).^2,'all')), data.syn1000(meta.csPCa_status==1),data.syn500(meta.csPCa_status==1),'UniformOutput', true);
%     
%     [~, p3] = ttest2(ben3, can3);
% end


end

