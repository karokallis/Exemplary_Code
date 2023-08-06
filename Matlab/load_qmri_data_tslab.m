function [outfile, data] = load_qmri_data_tslab(metadata, varargin)

%% Help for load_qmri_data.m

% Created: 03-24-22 by Tyler Seibert
% Adapted by Karoline Kallis

% Purpose: Load data for patient-level prostate RSI analysis, etc.
%
% Required Arguments:
%   subjects     =   excel list with patient IDs which correspond to folder
%   name -- toDo: might have to be updated based on redcap excel list

%%toDO: Select MRI according to data in Excel list

%% Process arguments:
defaults = {
    'organdir','/space/bil-syn01/1/cmig_bil/RSIData/Prostate/';
    'institutionname','UCSD';
    'datasetname','FPRSIv2';
    'studyIDprefix','PDS';
    'seriesname','Series_*_RSI_FOCUS_PROSTATE';
    'mapname','RSI_C_vol_NC';
    'segname', 'prostate_contour_DWI_space';
    'bladder_scalar_name','urine_norm_scalar';
    'mb0_name','mb0_scalar';
    'adcname','conventional_ADC_map'
    'dwiname','DWI_vol_NC_averaged'
    'dicomname','dcminfo'
    'bvaluesName','bvals'
    'urineName','urine_mask'
    'lesionname','*ROIs'
    'forceflag',0;
    'RepMeasurement',0; % if 1 = repeated measuremtn exists
    'outdir','';
    'outfile','';
    };
% procecdss variable inputs:
parse_varargin_auto(varargin,defaults)
%parse_varargin(varargin,defaults)
if isempty(outfile)
%     outfile = sprintf('%s%s_%s_%s_%s.mat',outdir,institutionname, datasetname,seriesname,mapname);
%     if contains(seriesname,'TEmin')
        outfile = sprintf('%s%s_%s_%s_%s_r%i.mat',outdir,institutionname,datasetname,seriesname,mapname,RepMeasurement);
%     end
end


%% Check if output exists:
if exist(outfile,'file') && forceflag==0
    fprintf('%s: output file already exists. Just returning file name: \n\t%s\n',mfilename,outfile);
    data = load(outfile);
    return
end

%% Add dependencies:
%daledir='/home/tylers/code/Dale_Matlab/matlab';
daledir='/space/bil-syn01/1/cmig_bil/code-repo/Tools/3_MghUtils';
pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(daledir, pathCell))
    fprintf('%s: Adding dependencies to path\n',mfilename);
    
end

%% initialize:

% nsubj = numel(subjects);
% outmaps = cell(1,numel(subjects));
% outseg = cell(1,numel(subjects));

outmaps = {};
outmaps_normalized = {};
outseg = {};
lesionseg = {};
incSubjects = {};
lname = {};

subjects = metadata.Study_ID;
dateList = datetime(metadata.MRI_Date, 'InputFormat','dd-MM-yyyy','Format','yyyyMMdd');
csPCa = metadata.csPCA;
gleason = metadata.GGG;
% hv_C = linspace(0,0.1,1001);
%% Get data:
k = 1;
for si=1:numel(subjects)
%for si=1:3
    
    map = {};
    rsi_TE = {};
    dwi={};
    seg ={};
    adc = {};
    mb0 = {};
    bladder_scalar = {};
    urine = {};
    subjects(si)
    substr = sprintf('%s%03.0f',char(subjects(si)));
    
    %substr = sprintf('%s%03.0f',studyIDprefix,char(subjects(si)));
    %     tmpstr = sprintf('%s/%s/%s/proc/%s/20*',organdir,institutionname,datasetname,substr);
    %     subjdir = get_ls(tmpstr,'-d',1);
    tmpstr = sprintf('%s/%s/PDS_MRI/proc/%s/%s/',organdir,institutionname,substr,dateList(si));
    subjdir = get_ls(tmpstr,'-d',1);
    %subjdir = {tmpstr};
    if isempty(subjdir)
        fprintf('%s: No data found for subject %s.\n\tLooked in %s\n',mfilename,substr,tmpstr);
    elseif length(subjdir)==1
        subjdir1 = subjdir{1};
        procdirstr = sprintf('%s/%s',subjdir1,seriesname);
        procdir = get_ls(procdirstr,'-d',1);
        
        if isempty(procdir)
            procdir1 = '';
            fprintf('%s: Could not find data for subject %s.\n\tLooked for %s\n',mfilename,substr,procdirstr);
        elseif length(procdir)==1
            procdir1 = procdir{1};
        elseif length(procdir)>1
            seriesNum = {};
            for l =1:length(procdir)
                tmpstr = split(procdir{l},'/');
                tmpstr = split(tmpstr{end},'_');
                seriesNum{l} = str2double(tmpstr{2});
           end
            if RepMeasurement==0 % might have to sort it by number that would actully be good
                [~, idx] = min(cell2mat(seriesNum));
                procdir1 = procdir{idx};
            else
                [~, idx] = max(cell2mat(seriesNum));
                procdir1 = procdir{idx}; % assuming if repeated, was because first was not ideal. May not be true.
            end
        end
        
        %% load outmaps --> store as outmaps{si}
        
        % load the RSI C maps:
        fname_mgz = sprintf('%s/%s.mgz',procdir1,mapname);
        dwiname_mgz = sprintf('%s/%s.mgz',procdir1,dwiname);
        sname_mgz = sprintf('%s/%s.mgz',procdir1,segname);
        lname_mat = sprintf('%s/%s.mat',procdir1,lesionname);
        lname = get_ls(lname_mat,'-d',1);
        
        adc_mgz = sprintf('%s/%s.mgz',procdir1,adcname);
        urine_mgz = sprintf('%s/%s.mgz',procdir1,urineName);
        dicom_mat = sprintf('%s/%s.mat',procdir1,dicomname);
        bvalues_mat = sprintf('%s/%s.mat',procdir1,bvaluesName);

        mb0_mat = sprintf('%s/%s.mat',procdir1,mb0_name);
        scalar_mat = sprintf('%s/%s.mat',procdir1,bladder_scalar_name);
                
        if ~exist(fname_mgz,'file') && ~exist(sname_mgz,'file')
            fprintf('%s: Could not find data or segmentation for subject %s. \n\tLooked for %s\n',mfilename,substr,fname_mgz);
        else
            map = QD_load_mgh(fname_mgz);
            dwi = QD_load_mgh(dwiname_mgz);
            seg = QD_load_mgh(sname_mgz);
            urine = QD_load_mgh(urine_mgz);
            adc = QD_load_mgh(adc_mgz);
            mb0 = load(mb0_mat);
            dicominfo = load(dicom_mat).dcminfo;
            TE = load(dicom_mat).dcminfo.EchoTime;
            bvalues = unique(load(bvalues_mat).bvals);
            bladder_scalar = load(scalar_mat);
            volume = nnz(seg)*prod([load(dicom_mat).dcminfo.PixelSpacing' load(dicom_mat).dcminfo.SliceThickness])/1000; %cc
            if contains(mapname,'RSI')
                map_normbladder = 100*map./bladder_scalar.k;
                rsi_TE = calculateRSI_TE(dwi,bvalues, TE, csPCa(si));
                %rsi_TE = calculateRSI_TE(dwi,bvalues, TE, csPCa(si))./bladder_scalar.k;
%                 vol_tmp = map_mb0_normbladder(:,:,:,1); % histogram only for C1 map
%                 [hc, ~] = hist(vol_tmp(seg>0.5),hv_C); % hv_C random binning // hc are counts // make histogram of voxel
%                 hc = hc/sum(hc);
                
            end
        end

        if ~isempty(lname)
            if ~exist(lname{1},'file')
                fprintf('%s: Could not find lesion contour for subject %s. \n\tLooked for %s\n',mfilename,substr,lesionname);
                lesion = [];
            else
                lesion = load(lname{1});
            end
        else
            lesion = [];
            
        end
        
        if(~isempty(map))
            dicomInfos{k} = dicominfo;
            outmaps{k} = map;
            outdwi{k} = dwi;
%             outmaps_TE{k} = rsi_TE;
            if(contains(mapname, 'RSI'))
                outmaps_normalized{k} = map_normbladder;
                outmaps_TE{k} = rsi_TE;
            end
            outseg{k} = seg ;
            ADCmap{k} = adc;
            urineMask{k} = urine;
            mb0_scalar{k}=mb0.k;
            normbladder_scalar{k} = bladder_scalar.k;
            csPCa_status{k} = csPCa(si);
            GGG{k} = gleason(si);
            balues{k} = bvalues;
            TEs{k} = TE;
            ProstateVolume{k} = volume;
            if ~isempty(lesion)
                lesionseg{k} = lesion;
            end
            incSubjects{k} = subjects{si};
            k = k+1;
        end
        clear fname_mgz
    end
end

% if ~isempty(lname)
%     data = struct('subjects', incSubjects, 'outmaps', outmaps, 'outmaps_TE',outmaps_TE,'outmaps_normalized',outmaps_normalized,'outseg', outseg, 'lesionseg', lesionseg,'adc',ADCmap,'mb0_scalar',mb0_scalar,'normbladder_scalar',normbladder_scalar);
% else
%     data = struct('subjects', incSubjects, 'outmaps', outmaps, 'outmaps_TE',outmaps_TE,'outmaps_normalized',outmaps_normalized,'outseg', outseg,'adc',ADCmap,'mb0_scalar',mb0_scalar,'normbladder_scalar',normbladder_scalar);
% end

subjects = incSubjects;
outfile = sprintf('%s%s_%s_%s_%s_r%i.mat',outdir,institutionname,datasetname,seriesname,mapname,RepMeasurement);

if ~isempty(lesionseg)
    save(outfile,'outmaps','outdwi','outmaps_normalized','outmaps_TE','outseg','urineMask','lesionseg','subjects','dicomInfos','ADCmap','mb0_scalar','normbladder_scalar','csPCa_status','GGG','bvalues','TEs','ProstateVolume','mapname','-v7.3');
else
    save(outfile,'outmaps','outdwi','outmaps_normalized','outmaps_TE','outseg','urineMask','subjects','dicomInfos','ADCmap','mb0_scalar','normbladder_scalar','csPCa_status','GGG','bvalues','TEs','mapname','ProstateVolume','-v7.3');
end

data = load(outfile);
fprintf('%s: Saved %s\n',mfilename,outfile);
