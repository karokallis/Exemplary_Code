function compareStats(datafiles,class,method, figName)
%Function to compare main stats between different sequences in datafiles
% mean
% median
% max
% min
% percentile 98
% percentile 99
% percentile 95



for i=1:length(datafiles)
    switch method
        case 'normalized'
            tmpcellMAX = cellfun(@(x,y) max(x(y>0.5),[],'all'), datafiles(i).data.outmaps_normalized(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            %tmpcellMIN = cellfun(@(x,y) min(nonzeros(x(y>0.5)),[],'all'), datafiles(i).data.outmaps_normalized(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellMEAN = cellfun(@(x,y) mean(x(y>0.5),'all'), datafiles(i).data.outmaps_normalized(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellMEDIAN = cellfun(@(x,y) median(x(y>0.5),'all'), datafiles(i).data.outmaps_normalized(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellPER95 = cellfun(@(x,y) prctile(x(y>0.5),95), datafiles(i).data.outmaps_normalized(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellPER98 = cellfun(@(x,y) prctile(x(y>0.5),98), datafiles(i).data.outmaps_normalized(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellPER99 = cellfun(@(x,y) prctile(x(y>0.5),99), datafiles(i).data.outmaps_normalized(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
        case 'TE'
            tmpcellMAX = cellfun(@(x,y) max(x(y>0.5),[],'all'), datafiles(i).data.outmaps_TE(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            %tmpcellMIN = cellfun(@(x,y) min(nonzeros(x(y>0.5)),[],'all'), datafiles(i).data.outmaps_normalized(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellMEAN = cellfun(@(x,y) mean(x(y>0.5),'all'), datafiles(i).data.outmaps_TE(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellMEDIAN = cellfun(@(x,y) median(x(y>0.5),'all'), datafiles(i).data.outmaps_TE(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellPER95 = cellfun(@(x,y) prctile(x(y>0.5),95), datafiles(i).data.outmaps_TE(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellPER98 = cellfun(@(x,y) prctile(x(y>0.5),98), datafiles(i).data.outmaps_TE(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellPER99 = cellfun(@(x,y) prctile(x(y>0.5),99), datafiles(i).data.outmaps_TE(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
        otherwise
            tmpcellMAX = cellfun(@(x,y) max(x(y>0.5),[],'all'), datafiles(i).data.outmaps(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            %tmpcellMIN = cellfun(@(x,y) min(nonzeros(x(y>0.5)),[],'all'), datafiles(i).data.outmaps_normalized(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellMEAN = cellfun(@(x,y) mean(x(y>0.5),'all'), datafiles(i).data.outmaps(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellMEDIAN = cellfun(@(x,y) median(x(y>0.5),'all'), datafiles(i).data.outmaps(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellPER95 = cellfun(@(x,y) prctile(x(y>0.5),95), datafiles(i).data.outmaps(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellPER98 = cellfun(@(x,y) prctile(x(y>0.5),98), datafiles(i).data.outmaps(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
            tmpcellPER99 = cellfun(@(x,y) prctile(x(y>0.5),99), datafiles(i).data.outmaps(:,:,:,1), datafiles(i).data.outseg, 'UniformOutput', false)';
  
    end
    
    switch class
        case 'cancer'
            idx{i} = cell2mat(datafiles(i).data.csPCa_status)==1;
        case 'benign'
            idx{i}= cell2mat(datafiles(i).data.csPCa_status)==0;
        otherwise
            tmp = cell2mat(tmpcellMAX);
            idx{i} = ones(1,length(tmp));
    end
    
    maxC1Prostate{i} = cell2mat(tmpcellMAX(idx{i}==1));
    %minC1Prostate{i} = cell2mat(tmpcellMIN(idx{i}==1));
    meanC1Prostate{i} = cell2mat(tmpcellMEAN(idx{i}==1));
    medianC1Prostate{i} = cell2mat(tmpcellMEDIAN(idx{i}==1));
    per95C1Prostate{i} = cell2mat(tmpcellPER95(idx{i}==1));
    per98C1Prostate{i} = cell2mat(tmpcellPER98(idx{i}==1));
    per99C1Prostate{i} = cell2mat(tmpcellPER99(idx{i}==1));
    
    %maxC1Prostate{i}  = rmoutliers(maxC1Prostate{i} , 'percentile',[0 99]);
end


%%

colors = [0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840;0 0 0.94];
grp = [repmat(["RSI FOCUS"], 1, length(maxC1Prostate{1})) repmat(["TE90"], 1, length(maxC1Prostate{2}))...
    repmat(["TEmin1"], 1, length(maxC1Prostate{3})) repmat(["TEmin2"], 1, length(maxC1Prostate{4}))];

fig = figure;
subplot(2,3,1)
C = [maxC1Prostate{1}; maxC1Prostate{2} ;maxC1Prostate{3} ;maxC1Prostate{4}];
boxplot(C,grp)
ylabel('Maximum SI')
%ylim([0,100])

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

% subplot(3,3,2)
% C = [minC1Prostate{1}; minC1Prostate{2} ;minC1Prostate{3} ;minC1Prostate{4}];
% grp = [repmat(["RSE"], 1, length(minC1Prostate{1})) repmat(["TE90"], 1, length(minC1Prostate{2}))...
%     repmat(["TEmin1"], 1, length(minC1Prostate{3})) repmat(["TEmin2"], 1, length(minC1Prostate{4}))];
% boxplot(C,grp)
% ylabel('Minimum SI')

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

subplot(2,3,2)
C = [meanC1Prostate{1}; meanC1Prostate{2} ;meanC1Prostate{3} ;meanC1Prostate{4}];
grp = [repmat(["RSI FOCUS"], 1, length(meanC1Prostate{1})) repmat(["TE90"], 1, length(meanC1Prostate{2}))...
    repmat(["TEmin1"], 1, length(meanC1Prostate{3})) repmat(["TEmin2"], 1, length(meanC1Prostate{4}))];

boxplot(C,grp)
%ylim([0,10])
ylabel('Mean SI')

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

subplot(2,3,3)

grp = [repmat(["RSI FOCUS"], 1, length(medianC1Prostate{1})) repmat(["TE90"], 1, length(medianC1Prostate{2}))...
    repmat(["TEmin1"], 1, length(medianC1Prostate{3})) repmat(["TEmin2"], 1, length(medianC1Prostate{4}))];
C = [medianC1Prostate{1}; medianC1Prostate{2} ;medianC1Prostate{3} ;medianC1Prostate{4}];
boxplot(C,grp)
ylabel('Median SI')
%ylim([0,10])

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

subplot(2,3,4)

grp = [repmat(["RSI FOCUS"], 1, length(per95C1Prostate{1})) repmat(["TE90"], 1, length(per95C1Prostate{2}))...
    repmat(["TEmin1"], 1, length(per95C1Prostate{3})) repmat(["TEmin2"], 1, length(per95C1Prostate{4}))];
C = [per95C1Prostate{1}; per95C1Prostate{2} ;per95C1Prostate{3} ;per95C1Prostate{4}];
boxplot(C,grp)
ylabel('95th Percentile SI')
%ylim([0,100])

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

subplot(2,3,5)
grp = [repmat(["RSI FOCUS"], 1, length(per98C1Prostate{1})) repmat(["TE90"], 1, length(per98C1Prostate{2}))...
    repmat(["TEmin1"], 1, length(per98C1Prostate{3})) repmat(["TEmin2"], 1, length(per98C1Prostate{4}))];
C = [per98C1Prostate{1}; per98C1Prostate{2} ;per98C1Prostate{3} ;per98C1Prostate{4}];
boxplot(C,grp)
ylabel('98th Percentile SI')
%ylim([0,100])

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end


subplot(2,3,6)

grp = [repmat(["RSI FOCUS"], 1, length(per99C1Prostate{1})) repmat(["TE90"], 1, length(per99C1Prostate{2}))...
    repmat(["TEmin1"], 1, length(per99C1Prostate{3})) repmat(["TEmin2"], 1, length(per99C1Prostate{4}))];
C = [per99C1Prostate{1}; per99C1Prostate{2} ;per99C1Prostate{3} ;per99C1Prostate{4}];
boxplot(C,grp)
ylabel('99th Percentile SI')
%ylim([0,100])


h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

sgtitle(figName)
%print(fig, ['images/Boxi_',figName],'-dpng','-r300');
end

