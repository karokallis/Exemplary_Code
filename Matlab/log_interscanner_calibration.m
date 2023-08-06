%% Script to run and analyse calibration between different series
%% toDo:
% make sure TEmin is sorted correctly
close all; clear all; clc;

% Created 03-21-23 Karoline Kallis
addpath(genpath('utils/'))
addpath(genpath('plots/'))
addpath(genpath('data/'))
addpath(genpath('../../synBvalues/matlab/ROC/'))

%addpath(genpath('/space/bil-syn01/1/cmig_bil/code-repo/Tools/utils/'))
addpath(genpath('/space/bil-syn01/1/cmig_bil/code-repo'))
% Remove this stupid thing that has a custom histogram function that has the same name as the matlab built-in function
rmpath(genpath('/space/bil-syn01/1/cmig_bil/code-repo/Tools/10_MorphometryToolbox/'))
%addpath(genpath('/home/dale/matlab/showVol_amd'));


%% NGP

% staticsubsfile = '/home/kkallis/Calibration/data/ngp/NGP_ptlevel_metadata.mat'; 
% tmpsublist = load(staticsubsfile).ptlevel;
% subjects = tmpsublist.Study_ID; % should replace prior subjects listtmp
try
    metadata = load('/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/NGP/NGP_ptlevel_metadata.mat').ptlevel;
catch
    disp('No saved metadata file found! Reading in metadata from raw file!')
    staticsubsfile = '/home/kkallis/RedCap/data/NGP_meta_26-Jul-2023.xlsx';
    tmpsublist = readtable(staticsubsfile, 'VariableNamingRule','preserve');
    metadata = load_meta_data_tslab(tmpsublist,'datasetname','NGP','studyIDprefix','PDS','outdir','/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/NGP/'); % PDS =KOP
end
% subjects = metadata.Study_ID;
% MRI_dates = datetime(metadata.MRI_Date, 'InputFormat','dd-MM-yyyy','Format','yyyyMMdd');
% csPCa = metadata.csPCa;
%% Load NGP series
datasetname = 'NGP';
studyIDprefix = 'PDS';
seriesname = '*RSI_FOCUS_PROSTATE';
outdir = '/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/NGP/';

%%
% RSI C maps, not normalized by urine
datafiles(1).seriesname = 'RSI_FOCUS_PROSTATE';
datafiles(1).dataset = datasetname;
[datafiles(1).fname, datafiles(1).data] = load_qmri_data_tslab(metadata, 'datasetname',datasetname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);


%% TE90
seriesname = '*FOCUS_TE90_b0_50_800_1500_3000';
% RSI C maps, not normalized by urine
datafiles(2).seriesname = 'Focus_TE90';
datafiles(2).dataset = datasetname;
[datafiles(2).fname, datafiles(2).data] = load_qmri_data_tslab(metadata,'datasetname',datasetname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);


%% TEmin!!!
seriesname = '*FOCUS_TEmin_b0_50_800_1500_3000';
% RSI C maps, not normalized by urine
datafiles(3).seriesname = 'Focus_TEmin';
datafiles(3).dataset = datasetname;
[datafiles(3).fname, datafiles(3).data] = load_qmri_data_tslab(metadata,'datasetname',datasetname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

%% TEmin!!! 2nd TEmin
seriesname = '*FOCUS_TEmin_b0_50_800_1500_3000';
% RSI C maps, not normalized by urine
datafiles(4).seriesname = 'Focus_TEmin';
datafiles(4).dataset = datasetname;
[datafiles(4).fname, datafiles(4).data] = load_qmri_data_tslab(metadata,'datasetname',datasetname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname, 'RepMeasurement',1);

%% Ensure that only cases are used which are in both datasets to be compared!
datafiles(2).data.outseg(86) = datafiles(3).data.outseg(86);
% something is up with the contour!!! 

%% Estimate scaling factor using linear Regression
[TEmin, TE90, TEmin_rep]= preprocessingDatasets(datafiles(3),datafiles(2),datafiles(4), 'benign');
%[~, TEmin_rep]= preprocessingDatasets(datafiles(3),datafiles(4), 'benign'); % Make sure that all cases are in both datasets
% Make sure that all cases are in both datasets
%[f, f_inter, R2 ,fall, ~] = estimateScalingFactor(TEmin.data,TE90.data,5,10,false);
f = load('data/scalingFactor_v2.mat');

%mb0_ref = mean([cell2mat(TEmin.data.mb0_scalar)',cell2mat(TEmin_rep.data.mb0_scalar)'],2);

MAD= compareSeries_v1(TE90, TEmin, TEmin_rep, f.sfvec, true,false);
mb0Factor_LR = compareRSIrs(TE90, TEmin, TEmin_rep, f.sfvec,'LinearRegression', true,false);
mb0Factor_T2 = compareRSIrs(TE90, TEmin, TEmin_rep, f.sfvec,'T2estimation', true,false);

%MAD_TEmin2_TEmin90= compareSeriesDWI(Ref, Comp, [1,1,1,1], true);
[TEmin, TE90, TEmin_rep]= preprocessingDatasets(datafiles(3),datafiles(2),datafiles(4), 'all');
MAD= compareSeries_v1(TE90, TEmin, TEmin_rep, f.sfvec, false,true);
% compareRSIrs(TE90, TEmin, TEmin_rep, f.sfvec,'input',mb0Factor_LR, true,true);
compareRSIrs(TE90, TEmin, TEmin_rep, f.sfvec,'T2estimation', false,true);






