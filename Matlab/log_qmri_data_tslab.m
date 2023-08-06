%% Log for patient-level analysis of FPRSI data
close all; clear all; clc;

% Created 03-09-22
% Adapted 04-13-22 Karoline Kallis
addpath(genpath('utils/'))
addpath(genpath('/space/bil-syn01/1/cmig_bil/code-repo/Tools/utils/'))
addpath(genpath('/space/bil-syn01/1/cmig_bil/code-repo'))
%% fprsi pre-processing data

% get subjects:
% % % % Still RedCap spreadsheet -- in case data is loaded from .xlsx
staticsubsfile = '/home/kkallis/RedCap/data/FPRSI_meta_09-May-2023.xlsx';
%staticsubsfile = '/home/kkallis/RedCap/data/FPRSI_meta_15-Feb-2023_n150.xlsx';
tmpsublist = readtable(staticsubsfile, 'VariableNamingRule','preserve');
metadata = load_meta_data_tslab(tmpsublist,'datasetname','FPRSI','studyIDprefix','PDS','outdir','/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/FPRSI/'); % PDS =KOP

%outfile = ('/home/kkallis/Calibration/data/fprsi/FPRSI_ptlevel_metadata_ptlevel_n151_pathCumulative.mat');
% staticsubsfile = '/home/kkallis/Calibration/data/DeIDed_Sheets/FPRSI_DeIDed.xlsx'; % fprsi151
% tmpsublist = readtable(staticsubsfile, 'VariableNamingRule','preserve');
%subjects = metadata.Study_ID; % should replace prior subjects listtmp
%dateList = datetime(metadata.MRI_Date, 'InputFormat','dd-MM-yyyy','Format','yyyyMMdd');
% metadata = load(outfile);
% subjects = metadata.metasubjectIDs; 
%metadata =  (tmpsublist,'datasetname','FPRSI');

datasetname = 'FPRSIv2';
studyIDprefix = 'PDS';
seriesname = 'Series_*__RSI_FOCUS_PROSTATE';
outdir = '/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/FPRSI/';
%%
% DWI pre noise correction
datafiles(1).dataset = datasetname;
datafiles(1).mapname = 'DWI_vol_preprocessed_averaged';
[datafiles(1).fname, datafiles(1).data] = load_qmri_data_tslab(metadata,'datasetname',datasetname,'mapname',datafiles(1).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);


% DWI post-noise correction
datafiles(2).dataset = datasetname;
datafiles(2).mapname = 'DWI_vol_NC_averaged';
[datafiles(2).fname, datafiles(2).data] = load_qmri_data_tslab(metadata,'datasetname',datasetname,'mapname',datafiles(2).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);
%%
% RSI C maps, not normalized by urine
datafiles(3).mapname = 'RSI_C_vol_NC';
datafiles(3).dataset = datasetname;
[datafiles(3).fname, datafiles(3).data] = load_qmri_data_tslab(metadata,'datasetname',datasetname,'mapname',datafiles(3).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);
%%
% % RSI C maps after normalization by urine
% datafiles(4).mapname = 'RSI_C_vol_normBladder';
% datafiles(4).dataset = datasetname;
% [datafiles(4).fname, datafiles(4).data]  = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(4).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% RSI C maps after normalization by urine and without noise correction
datafiles(4).mapname = 'RSI_C_vol_noNC';
datafiles(4).dataset = datasetname;
[datafiles(4).fname, datafiles(4).data]  = load_qmri_data_tslab(metadata,'datasetname',datasetname,'mapname',datafiles(4).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

%% DREAM
% Dream pre-processing data

% staticsubsfile = '/home/kkallis/Calibration/data/dream/MEPR_ptlevel_metadata.mat'; 
% tmpsublist = load(staticsubsfile).ptlevel;
% subjects = tmpsublist.StudyID; % should replace prior subjects listtmp


staticsubsfile = '/home/kkallis/RedCap/data/DREAM_meta_30-Sep-2022.xlsx';
tmpsublist = readtable(staticsubsfile, 'VariableNamingRule','preserve');
metadata = load_meta_data_tslab(tmpsublist,'datasetname','MEPR','studyIDprefix','PDS','outdir','/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/Dream/'); % PDS =KOP
dateList = datetime(metadata.MRI_Date, 'InputFormat','dd-MM-yyyy','Format','yyyyMMdd');

% % Still Christines spreadsheet -- in case data is loaded from .xlsx
% staticsubsfile = '/home/kkallis/Calibration/data/DeIDed_Sheets/Dream_Meta_v1.xlsx';
% tmpsublist = readtable(staticsubsfile, 'VariableNamingRule','preserve');
% metadata = load_meta_data_tslab(tmpsublist,'datasetname','MEPR','outdir','/home/kkallis/Calibration/data/dream'); % PDS =KOP
subjects = metadata.Study_ID; % should replace prior subjects listtmp
%%
datasetname = 'Dreamv2';
studyIDprefix = 'PDS';
seriesname = 'Series_*__MMIL_MULTISHELL_DIFFUSION';
outdir = '/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/Dream/';

% DWI pre noise correction
datafiles(5).dataset = datasetname;
datafiles(5).mapname = 'DWI_vol_preprocessed_averaged';
[datafiles(5).fname, datafiles(5).data] = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(5).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% DWI post-noise correction
datafiles(6).dataset = datasetname;
datafiles(6).mapname = 'DWI_vol_NC_averaged';
[datafiles(6).fname, datafiles(6).data] = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(6).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% RSI C maps, not normalized by urine
datafiles(7).mapname = 'RSI_C_vol_NC';
datafiles(7).dataset = datasetname;
[datafiles(7).fname, datafiles(7).data] = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(7).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% 
% % RSI C maps after normalization by urine
% datafiles(9).mapname = 'RSI_C_vol_normBladder';
% datafiles(9).dataset = datasetname;
% [datafiles(9).fname, datafiles(9).data]  = load_qmri_data_tslab(subjects,'datasetname',datasetname,'mapname',datafiles(9).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);


% RSI C maps after normalization by urine and without noise correction
datafiles(8).mapname = 'RSI_C_vol_noNC';
datafiles(8).dataset = datasetname;
[datafiles(8).fname, datafiles(8).data]  = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(8).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

%% load data for these subjects: TE_80
seriesname = '*MMIL_*_TE80';


% DWI pre noise correction
datafiles(9).dataset = datasetname;
datafiles(9).mapname = 'DWI_vol_preprocessed_averaged';
[datafiles(9).fname, datafiles(9).data] = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(9).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% DWI post-noise correction
datafiles(10).dataset = datasetname;
datafiles(10).mapname = 'DWI_vol_NC_averaged';
[datafiles(10).fname, datafiles(10).data] = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(10).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% RSI C maps, not normalized by urine
datafiles(11).mapname = 'RSI_C_vol_NC';
datafiles(11).dataset = datasetname;
[datafiles(11).fname, datafiles(11).data] = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(11).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% 
% % RSI C maps after normalization by urine
% datafiles(9).mapname = 'RSI_C_vol_normBladder';
% datafiles(9).dataset = datasetname;
% [datafiles(9).fname, datafiles(9).data]  = load_qmri_data_tslab(subjects,'datasetname',datasetname,'mapname',datafiles(9).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);


% RSI C maps after normalization by urine and without noise correction
datafiles(12).mapname = 'RSI_C_vol_noNC';
datafiles(12).dataset = datasetname;
[datafiles(12).fname, datafiles(12).data]= load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(12).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);


%% load RSI data for these subjects: TE_100

seriesname = '*MMIL_*_TE100';
% DWI pre noise correction
% DWI pre noise correction
datafiles(13).dataset = datasetname;
datafiles(13).mapname = 'DWI_vol_preprocessed_averaged';
[datafiles(13).fname, datafiles(13).data] = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(13).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% DWI post-noise correction
datafiles(14).dataset = datasetname;
datafiles(14).mapname = 'DWI_vol_NC_averaged';
[datafiles(14).fname, datafiles(14).data] = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(14).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% RSI C maps, not normalized by urine
datafiles(15).mapname = 'RSI_C_vol_NC';
datafiles(15).dataset = datasetname;
[datafiles(15).fname, datafiles(15).data] = load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(15).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% 
% % RSI C maps after normalization by urine
% datafiles(9).mapname = 'RSI_C_vol_normBladder';
% datafiles(9).dataset = datasetname;
% [datafiles(9).fname, datafiles(9).data]  = load_qmri_data_tslab(subjects,'datasetname',datasetname,'mapname',datafiles(9).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);


% RSI C maps after normalization by urine and without noise correction
datafiles(16).mapname = 'RSI_C_vol_noNC';
datafiles(16).dataset = datasetname;
[datafiles(16).fname, datafiles(16).data]= load_qmri_data_tslab(subjects,dateList,'datasetname',datasetname,'mapname',datafiles(16).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);


%% NGP

% staticsubsfile = '/home/kkallis/Calibration/data/ngp/NGP_ptlevel_metadata.mat'; 
% tmpsublist = load(staticsubsfile).ptlevel;
% subjects = tmpsublist.Study_ID; % should replace prior subjects listtmp

% % % % Still RedCap spreadsheet -- in case data is loaded from .xlsx
staticsubsfile = '/home/kkallis/RedCap/data/NGP_ProRSI_meta_01-Mar-2023.xlsx';
tmpsublist = readtable(staticsubsfile, 'VariableNamingRule','preserve');
metadata = load_meta_data_tslab(tmpsublist,'datasetname','NGP','studyIDprefix','PDS','outdir','/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/NGP/'); % PDS =KOP
subjects = metadata.Study_ID;


%% ProRSI
staticsubsfile = '/home/kkallis/RedCap/data/ProRSI_meta_02-Mar-2023.xlsx';
tmpsublist = readtable(staticsubsfile, 'VariableNamingRule','preserve');
metadata_ProRSI = load_meta_data_tslab(tmpsublist,'datasetname','ProRSI','studyIDprefix','PDS','outdir','/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/NGP/'); % PDS =KOP

% 
subjects_ProRSI = metadata_ProRSI.Study_ID;
% 
% idx = ~ismember(subjects_ProRSI,subjects_NGP);
% %%
% 
% metadata_ProRSI(:,[2]) = [];% delete ProRSI ID
% metadata.Lesion_zone = [];
% metadata_ProRSI.Lesion_zone = [];
% 
% metadata.Lesion_Size = [];
% metadata_ProRSI.Lesion_Size = [];
% 
% metadata.psa2mri = [];
% metadata_ProRSI.psa2mri = [];
% 
% ptlevel = [metadata; metadata_ProRSI];
% 
% subjects = [subjects_NGP; subjects_ProRSI(idx)];
% MRI_dates = [metadata.MRI_Date; metadata_ProRSI.MRI_Date(idx)];
MRI_dates = datetime(metadata.MRI_Date, 'InputFormat','dd-MM-yyyy','Format','yyyyMMdd');
%%
datasetname = 'NGP';
studyIDprefix = 'PDS';
seriesname = '*RSI_FOCUS_PROSTATE';
outdir = '/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/NGP/';

%%
% DWI pre noise correction
datafiles(17).dataset = datasetname;
datafiles(17).mapname = 'DWI_vol_preprocessed_averaged';
[datafiles(17).fname, datafiles(17).data] = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(17).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname,'lesionname','');
%%
%
% DWI post-noise correction
datafiles(18).dataset = datasetname;
datafiles(18).mapname = 'DWI_vol_NC_averaged';
[datafiles(18).fname, datafiles(18).data] = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(18).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% RSI C maps, not normalized by urine
datafiles(19).mapname = 'RSI_C_vol_NC';
datafiles(19).dataset = datasetname;
[datafiles(19).fname, datafiles(19).data] = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(19).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% % RSI C maps after normalization by urine
% datafiles(24).mapname = 'RSI_C_vol_normBladder';
% datafiles(24).dataset = datasetname;
% [datafiles(24).fname, datafiles(24).data]  = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(24).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% RSI C maps after normalization by urine and without noise correction
datafiles(20).mapname = 'RSI_C_vol_noNC';
datafiles(20).dataset = datasetname;
[datafiles(20).fname, datafiles(20).data]  = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(20).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);


%% TE90
seriesname = '*FOCUS_TE90_b0_50_800_1500_3000';


datafiles(21).dataset = datasetname;
datafiles(21).mapname = 'DWI_vol_NC_averaged';
[datafiles(21).fname, datafiles(21).data] = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(21).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);
%
% RSI C maps, not normalized by urine
datafiles(22).mapname = 'RSI_C_vol_NC';
datafiles(22).dataset = datasetname;
[datafiles(22).fname, datafiles(22).data] = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(22).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% % RSI C maps after normalization by urine
% datafiles(24).mapname = 'RSI_C_vol_normBladder';
% datafiles(24).dataset = datasetname;
% [datafiles(24).fname, datafiles(24).data]  = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(24).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% RSI C maps after normalization by urine and without noise correction
datafiles(23).mapname = 'RSI_C_vol_noNC';
datafiles(23).dataset = datasetname;
[datafiles(23).fname, datafiles(23).data]  = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(23).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% DWI pre noise correction
datafiles(24).dataset = datasetname;
datafiles(24).mapname = 'DWI_vol_preprocessed_averaged';
[datafiles(24).fname, datafiles(24).data] = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(24).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname,'lesionname','');


%% TEmin!!!
seriesname = '*FOCUS_TEmin_b0_50_800_1500_3000';
%%
% DWI pre noise correction
datafiles(25).dataset = datasetname;
datafiles(25).mapname = 'DWI_vol_preprocessed_averaged';
[datafiles(25).fname, datafiles(25).data] = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(25).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname,'lesionname','');

%%
% DWI post-noise correction
datafiles(26).dataset = datasetname;
datafiles(26).mapname = 'DWI_vol_NC_averaged';
[datafiles(26).fname, datafiles(26).data] = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(26).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% RSI C maps, not normalized by urine
datafiles(27).mapname = 'RSI_C_vol_NC';
datafiles(27).dataset = datasetname;
[datafiles(27).fname, datafiles(27).data] = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(27).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% % RSI C maps after normalization by urine
% datafiles(24).mapname = 'RSI_C_vol_normBladder';
% datafiles(24).dataset = datasetname;
% [datafiles(24).fname, datafiles(24).data]  = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(24).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

% RSI C maps after normalization by urine and without noise correction
datafiles(28).mapname = 'RSI_C_vol_noNC';
datafiles(28).dataset = datasetname;
[datafiles(28).fname, datafiles(28).data]  = load_qmri_data_tslab(subjects,MRI_dates,'datasetname',datasetname,'mapname',datafiles(28).mapname,'studyIDprefix',studyIDprefix,'outdir',outdir,'seriesname',seriesname);

%% KOP
staticsubsfile = '/home/kkallis/RedCap/data/KOP_meta_25-May-2023.xlsx';
tmpsublist = readtable(staticsubsfile, 'VariableNamingRule','preserve');
metadata = load_meta_data_tslab(tmpsublist,'datasetname','KOP','studyIDprefix','PDS','outdir','/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/Calibration/KOP/'); % PDS =KOP
%subjects_NGP = metadata.Study_ID;





