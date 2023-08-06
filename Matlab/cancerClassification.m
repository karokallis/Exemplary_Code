function [canceridx,benignidx,includeix] = cancerClassification(ptLevel,dataset)
% Created by Karoline Kallis 7/19/2022
% Function to classify cancer vs no cancer!
includepath = {'1','2','3','4','5','Benign'}; % Find relevant patients
includeidx = ismember(ptLevel.path,includepath)|ismember(ptLevel.WorstPIRADS,1); 
ptLevel = ptLevel(includeidx,:);

% Sort patient data to relevenat patients
includeix = ismember(dataset.subjects,ptLevel.StudyID);
dataset.outmaps = dataset.outmaps(includeix);
dataset.outseg = dataset.outseg(includeix);
dataset.subjects = dataset.subjects(includeix);


includedix = ismember(ptLevel.StudyID, dataset.subjects);
ptLevel = ptLevel(includedix,:);

%% Create 2 groups of patients cancer no cancer
% no cancer = pirads =1 / Gleason <=2?
% cancer = Gleason >3
% path = summary of Gleason score!
%

cancerpath = {'2','3','4','5'};
benignpath = {'Benign','1'};
canceridx = ismember(ptLevel.path,cancerpath); 
benignidx = ismember(ptLevel.path,benignpath); 

disp('BENIGN');
summary(ptLevel(benignidx,{'WorstPIRADS','path'}))


disp('CANCER');
summary(ptLevel(canceridx,{'WorstPIRADS','path'}))
% cancer = getMapCell(dataset.outmaps(canceridx),1);
% benign = getMapCell(dataset.outmaps(benignidx),1);

end

