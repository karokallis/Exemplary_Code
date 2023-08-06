function [dataset1Updated,dataset2Updated, dataset3Updated] = preprocessingDatasets(dataset1,dataset2,dataset3,class)
% Created by Karoline Kallis 05/12/2023
%Function to make sure that all datasets exist in the same datafiles

if isempty(dataset3)
    switch class
        case 'benign'
            idxComp = ismember(dataset1.data.subjects, dataset2.data.subjects)& cell2mat(dataset1.data.csPCa_status)==0;
            idxReference = ismember(dataset2.data.subjects, dataset1.data.subjects)& cell2mat(dataset2.data.csPCa_status)==0;
        case 'cancer'
            idxComp = ismember(dataset1.data.subjects, dataset2.data.subjects)& cell2mat(dataset1.data.csPCa_status)==1;
            idxReference = ismember(dataset2.data.subjects, dataset1.data.subjects)& cell2mat(dataset2.data.csPCa_status)==1;
        otherwise
            idxComp = ismember(dataset1.data.subjects, dataset2.data.subjects);
            idxReference = ismember(dataset2.data.subjects, dataset1.data.subjects);
    end
else
    A1 = string(setdiff(dataset1.data.subjects, dataset2.data.subjects))';
    A2 = string(setdiff(dataset2.data.subjects, dataset1.data.subjects))';

    B1 = string(setdiff(dataset1.data.subjects, dataset3.data.subjects))';
    B2 = string(setdiff(dataset3.data.subjects, dataset1.data.subjects))';

    C1 = string(setdiff(dataset2.data.subjects, dataset3.data.subjects))';
    C2 = string(setdiff(dataset3.data.subjects, dataset2.data.subjects))';

    
    D = unique([A1;A2;B1;B2;C1;C2]);
    switch class
        
        case 'benign'
            idxComp = ~ismember(dataset1.data.subjects, D)& cell2mat(dataset1.data.csPCa_status)==0;
            idxReference = ~ismember(dataset2.data.subjects, D)& cell2mat(dataset2.data.csPCa_status)==0;
            idxRep = ~ismember(dataset3.data.subjects, D)& cell2mat(dataset3.data.csPCa_status)==0;
        case 'cancer'
            idxComp = ~ismember(dataset1.data.subjects, D)& cell2mat(dataset1.data.csPCa_status)==1;
            idxReference = ~ismember(dataset2.data.subjects, D)& cell2mat(dataset2.data.csPCa_status)==1;
            idxRep = ~ismember(dataset3.data.subjects, D)& cell2mat(dataset3.data.csPCa_status)==1;
        otherwise
            idxComp = ~ismember(dataset1.data.subjects, D);
            idxReference = ~ismember(dataset2.data.subjects, D);
            idxRep = ~ismember(dataset3.data.subjects, D);
    end
end
dataset1Updated.data.ADCmap = dataset1.data.ADCmap(idxComp==1);
dataset1Updated.data.GGG = dataset1.data.GGG(idxComp==1);
dataset1Updated.data.bvalues = dataset1.data.bvalues;
dataset1Updated.data.ProstateVolume = dataset1.data.ProstateVolume(idxComp==1);
dataset1Updated.data.TEs = dataset1.data.TEs(idxComp==1);
dataset1Updated.data.csPCa_status = dataset1.data.csPCa_status(idxComp==1);
dataset1Updated.data.mb0_scalar = dataset1.data.mb0_scalar(idxComp==1);
dataset1Updated.data.normbladder_scalar = dataset1.data.normbladder_scalar(idxComp==1);
dataset1Updated.data.outdwi = dataset1.data.outdwi(idxComp==1);
dataset1Updated.data.outmaps = dataset1.data.outmaps(idxComp==1);
dataset1Updated.data.urineMask = dataset1.data.urineMask(idxComp==1);
dataset1Updated.data.outmaps_TE = dataset1.data.outmaps_TE(idxComp==1);
dataset1Updated.data.outmaps_normalized = dataset1.data.outmaps_normalized(idxComp==1);
dataset1Updated.data.outseg = dataset1.data.outseg(idxComp==1);
dataset1Updated.data.subjects = dataset1.data.subjects(idxComp==1);
dataset1Updated.data.dicomInfos = dataset1.data.dicomInfos(idxComp==1);

dataset2Updated.data.ADCmap = dataset2.data.ADCmap(idxReference==1);
dataset2Updated.data.GGG = dataset2.data.GGG(idxReference==1);
dataset2Updated.data.bvalues = dataset2.data.bvalues;
dataset2Updated.data.ProstateVolume = dataset2.data.ProstateVolume(idxReference==1);
dataset2Updated.data.TEs = dataset2.data.TEs(idxReference==1);
dataset2Updated.data.csPCa_status = dataset2.data.csPCa_status(idxReference==1);
dataset2Updated.data.mb0_scalar = dataset2.data.mb0_scalar(idxReference==1);
dataset2Updated.data.normbladder_scalar = dataset2.data.normbladder_scalar(idxReference==1);
dataset2Updated.data.outdwi = dataset2.data.outdwi(idxReference==1);
dataset2Updated.data.outmaps = dataset2.data.outmaps(idxReference==1);
dataset2Updated.data.outmaps_TE = dataset2.data.outmaps_TE(idxReference==1);
dataset2Updated.data.urineMask = dataset2.data.urineMask(idxReference==1);
dataset2Updated.data.outmaps_normalized = dataset2.data.outmaps_normalized(idxReference==1);
dataset2Updated.data.outseg = dataset2.data.outseg(idxReference==1);
dataset2Updated.data.subjects = dataset2.data.subjects(idxReference==1);
dataset2Updated.data.dicomInfos = dataset2.data.dicomInfos(idxReference==1);


dataset3Updated.data.ADCmap = dataset3.data.ADCmap(idxRep==1);
dataset3Updated.data.GGG = dataset3.data.GGG(idxRep==1);
dataset3Updated.data.bvalues = dataset3.data.bvalues;
dataset3Updated.data.ProstateVolume = dataset3.data.ProstateVolume(idxRep==1);
dataset3Updated.data.TEs = dataset3.data.TEs(idxRep==1);
dataset3Updated.data.csPCa_status = dataset3.data.csPCa_status(idxRep==1);
dataset3Updated.data.mb0_scalar = dataset3.data.mb0_scalar(idxRep==1);
dataset3Updated.data.normbladder_scalar = dataset3.data.normbladder_scalar(idxRep==1);
dataset3Updated.data.outdwi = dataset3.data.outdwi(idxRep==1);
dataset3Updated.data.outmaps = dataset3.data.outmaps(idxRep==1);
dataset3Updated.data.outmaps_TE = dataset3.data.outmaps_TE(idxRep==1);
dataset3Updated.data.urineMask = dataset3.data.urineMask(idxRep==1);
dataset3Updated.data.outmaps_normalized = dataset3.data.outmaps_normalized(idxRep==1);
dataset3Updated.data.outseg = dataset3.data.outseg(idxRep==1);
dataset3Updated.data.subjects = dataset3.data.subjects(idxRep==1);
dataset3Updated.data.dicomInfos = dataset3.data.dicomInfos(idxRep==1);

%DWI2 = dataset2Updated.data.outdwi;

%DWI1 = dataset1Updated.data.outdwi;
%DWI2 = cellfun(@(x,y) x./y, dataset2Updated.data.outdwi, dataset2Updated.data.normbladder_scalar, 'UniformOutput', false);

%DWI1 = cellfun(@(x,y) x./y, dataset1Updated.data.outdwi, dataset1Updated.data.normbladder_scalar, 'UniformOutput', false);
%% Make a sanity check that the correct subjects are included

for subji=1:length(dataset2Updated.data.subjects)
    if(~strcmp(dataset2Updated.data.subjects{subji},dataset1Updated.data.subjects{subji}))
        disp('Error subjects mismatch!')
        dataset2Updated.data.subjects{subji}
        dataset1Updated.data.subjects{subji}
    end
end

end

