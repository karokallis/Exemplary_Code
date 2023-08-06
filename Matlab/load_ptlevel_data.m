function outfile = load_ptlevel_data(varargin)

%% Help for qMRI_analysis.m

% Created: 03-25-22 by Tyler Seibert
% Adapted: 04-13-22 by Karoline Kallis
%   See Change Log at end of file for modifications
% 
% Purpose: 
%   -perform quantitative MRI analysis on prostate MRI data
% 
% Required Arguments:
%   -rawmetafile       File name for .mat file containing outmaps and subjects
% 

disp('');

%% Process arguments:

% if nargin < 1; error('%s: Must specify procdir\n',mfilename); end

defaults = {    
    'staticsubsfile','FPRSI_ptlevel_metadata_ptlevel_n151_pathCumulative.mat'; % fprsi151

    %'workdir','/space/syn02/1/data/MMILDB/tylers/RSI/fprsi';
    'workdir','/home/kkallis/Calibration/data/fprsi'; %% Might lead to some issues!!!
    'datadir',''; % will be named based on workdir
    'rawmetafile',''; % will be named based on datadir
    'contourdir','';  % will be named based on workdir
    'datasetname','FPRSI';
    'outdir','/home/kkallis/synBalvues/data';
    'outfile',''; 
        
    };
        
% process variable inputs:        
parse_varargin_auto(varargin,defaults)

if exist(staticsubsfile,'file'); static_subjects = 1; else static_subjects = 0; end % mat file doesnt exist
if isempty(datadir); datadir = sprintf('%s/datafiles',workdir); end
if isempty(outdir); outdir = sprintf('%s/%s',datadir,outstr); end
if isempty(rawmetafile); rawmetafile = sprintf('%s/%s_ptlevel_metadata.mat',datadir, datasetname); end
if isempty(contourdir); contourdir = sprintf('%s/mimcontours',workdir); end

% outdir:
if ~exist(outdir,'dir'); unxcmd(['mkdir ' outdir]); end

%% Check for dale code:
if ~exist('sfigure'); fprintf('adding dalematlab to path'); dalematlab; end

%% Get the data:

% metadata:
load(rawmetafile); % loads ptlevel
summary(ptlevel(:,{'WorstGleason','WorstPIRADS','path'}))
%summary(ptlevel(:,{'WorstGleason','WorstPIRADS'}))


if static_subjects
    % static subject list
    tmpsublist = readtable(staticsubsfile);
    subjects = tmpsublist.StudyID; % should replace prior subjects list
    includeidx = ismember(ptlevel.StudyID,subjects);
else
    includepath = {'1','2','3','4','5','Benign'}; % skip NoBiopsy and Other
    includepirads = {'1','2','3','4','5'}; % skip 0 (not sure what that even means
    includecontours = get_ls(contourdir); % contours available
    includeidx = ismember(ptlevel.path,includepath) ...
             & ismember(ptlevel.WorstPIRADS,includepirads) ...
             & ~ptlevel.TreatmentBeforeMRI ... % no treatment prior to MRI
             & abs(ptlevel.DaysMRIBiopsy)<=180 ... % biopsy within 180 days of MRI
             ;             %& ismember(ptlevel.StudyID,includecontours) ...
    subjects = ptlevel.subnum(includeidx);
end
path = ptlevel.path(includeidx);
pirads = ptlevel.WorstPIRADS(includeidx);
csPCa = {'2','3','4','5'};
csPCa_status = double(ismember(path,csPCa)); % convert to numeric because tms_plot_ROC.m and perfcurve.m do not accept logical 
summary(ptlevel(includeidx,{'path','WorstPIRADS'}))
sum(csPCa_status); % 2020.12.15: N=151. 83 csPCa, 39 GGG 1, 29 benign. 2020.12.11 (did not exclude prior Tx or DaysBiopsyMRI>180: 123 clinically sig cases, 55 GGG 1 cases, 41 benign 
%worstsectors = ptlevel.WorstSectors(includeidx);

% data by zone:
zone = ptlevel.Zone(includeidx);
tz_pirads = categorical(ptlevel.PR_WorstTZ(includeidx));
pz_pirads = categorical(ptlevel.PR_WorstPZ(includeidx));
subsetidx_tz = ismember(zone,{'T'}); % include if TZ lesion identified by PI-RADS
subsetidx_pz = ismember(zone,{'P'}); % include if PZ lesion identified by PI-RADS
% Note: Excluding those with PI-RADS lesions identified in both TZ and PZ
% because we do not have "worstGleason" broken down by TZ and PZ. We could
% have also included the benign/GGG1 biopsies, but we cannot be sure that
% the TZ was well sampled in those with PI-RADS 1 TZ these (definitely was 
% not in those that had MR-only biopsy, for example).
subsetidx_all = true(numel(subjects),1); % include all

% convert pirads to numeric (was categorical because we had 'NoBiopsy', etc.)
% Not necessary for the new DeID I got
% for si=1:numel(subjects)
%     all_piradsnum(si,1) = str2num(char(pirads(si)));
%     tz_piradsnum(si,1) = str2num(char(tz_pirads(si)));
%     pz_piradsnum(si,1) = str2num(char(pz_pirads(si)));
% end

metasubjects = subjects; % rename so it is different from the datafile subjects variable 

% save this patient-level metadata:
if isempty(outfile)
    tmpstr = sprintf('_ptlevel_n%d.mat',numel(subjects));
    outfile = regexprep(rawmetafile,'.mat',tmpstr)
end
if exist(outfile,'file')
    fprintf('%s: Output file already exists. \n\tNot re-writing %s\n',mfilename,outfile);
else
    save(outfile,'metasubjects','path','pirads','csPCa','csPCa_status','zone',...
        'tz_pirads','pz_pirads','subsetidx_tz','subsetidx_pz','subsetidx_all',...
        '-v7.3');%'all_piradsnum','tz_piradsnum','pz_piradsnum',...
        
    fprintf('%s: Saved %s\n',mfilename,outfile);
end

%% Change Log
% 03-25-22 --> Created by Tyler Seibert (TMS)
% 03-25-22 --> (TMS) Replaces beginning of qMRI_analysis.m. Run this first.
% 
% 
