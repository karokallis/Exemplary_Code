function metadata = load_meta_data_tslab(ptlevel,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Process arguments:
defaults = {
    'datasetname','FPRSI_v2';
    'studyIDprefix','PDS';
    'forceflag',0;
    'datadirpath','/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/PDS_MRI/proc/'
    'outdir','';
    'outfile',''; 
    };

% procecdss variable inputs:         
parse_varargin_auto(varargin,defaults)
%parse_varargin(varargin,defaults)

%% Reduce table to interesint inputs

% idx=contains(ptlevel.StudyID,datasetname);
% ptlevel = ptlevel(idx,:);
% ptlevel.path = convertGS_to_GGG(ptlevel.WorstGleason);
%% Check if output exists:
if exist(outfile,'file') && forceflag==0
    fprintf('%s: output file already exists. Just returning file name: \n\t%s\n',mfilename,outfile);
    return
end

if (isempty(outfile)||forcelflag ==1) 
    outfile = sprintf('%s/%s_%s.mat',outdir,datasetname,'ptlevel_metadata');
    ptlevel.path = convertGS_to_GGG(ptlevel.Biopsy_highGleason);
    ptlevel.hist = convertGS_to_GGG(ptlevel.prostectomy_Gleason);
    
    for i = 1:length(ptlevel.path)
        if(contains(string(ptlevel.path(i)),'NoBiopsy'))
            if(any(strcmp(ptlevel.Properties.VariableNames,'HighestPirads')))
                if ptlevel.HighestPirads(i) <=2
                    ptlevel.path(i) = 'pirads';
                elseif ~contains(string(ptlevel.hist(i)),'Other')
                    ptlevel.path(i) = ptlevel.hist(i);
                end
            else
                if ptlevel.Pirads(i) <=2
                    ptlevel.path(i) = 'pirads';
                elseif ~contains(string(ptlevel.hist(i)),'Other')
                    ptlevel.path(i) = ptlevel.hist(i);
                end
            end
            
        end
    end
    %% Add csPCA_status
    cancerpath = {'2','3','4','5'};
    %benignpath = {'Benign','1'};
    if(strcmp(datasetname,'ProRSI')) % should be updated! 
        ptlevel.csPCA = ones(length(ptlevel.path),1); % change once RedCap is updated!!!
    else
        ptlevel.csPCA = ismember(ptlevel.path,cancerpath)|ismember(ptlevel.hist,cancerpath); % 0 benign // 1 cancer
    end
    ptlevel.GGG = ptlevel.path;
    ptlevel.GGG(find(ptlevel.hist ~= 'Other')) = ptlevel.path(find(ptlevel.hist ~= 'Other'));
    
    ptlevel.GGG(ismember(ptlevel.GGG,{'Benign' 'Other','pirads'})==1) = '0';
    ptlevel.dataExist = zeros(length(ptlevel.csPCA),1);
    for su=1:height(ptlevel)
        id = num2str(ptlevel.Subject_ID(su));
        ptlevel.Study_ID(su) = strcat('PDS',convertCharsToStrings(id));
%         if(length(id)<4&&length(id)>2)
%             ptlevel.Study_ID(su) = strcat('PDS0',convertCharsToStrings(id));
%         elseif(length(id)<=2)
%             ptlevel.Study_ID(su) = strcat('PDS00',convertCharsToStrings(id));
%         else
%             ptlevel.Study_ID(su) = strcat('PDS',convertCharsToStrings(id));
%         end
        %% Check if folder exists to be processed
        MRI_dates = datetime(ptlevel.MRI_Date, 'InputFormat','dd-MM-yyyy','Format','yyyyMMdd');
        switch studyIDprefix
            case 'FPRSI'
                tmpstr = sprintf('%s%s/%s/',datadirpath,ptlevel.FPRSI_ID{su,1},MRI_dates(su));
                tmp = exist(tmpstr,'dir');
            case 'MEPR'
                tmpstr = sprintf('%s%s/%s/',datadirpath,ptlevel.DREAM_ID{su,1},MRI_dates(su));
                tmp = exist(tmpstr,'dir');
            otherwise
                tmpstr = sprintf('%s%s/%s/',datadirpath,ptlevel.Study_ID(su),MRI_dates(su));
                tmp = exist(tmpstr,'dir');
        end
        if tmp>0
            ptlevel.dataExist(su)=1;
            
            F = dir(tmpstr);
            fld = extractfield(F, 'name');
            
            ptlevel.imageSeries(su) = string(strjoin(fld(3:end),'; '));

        end
                
    end
%     prefix = repmat(studyIDprefix,height(ptlevel),1);
%     ptlevel.Study_ID = strcat(prefix,num2str(ptlevel.Subject_ID));
%     ptlevel.Study_ID = regexprep(ptlevel.Study_ID,' ','0');



end

save(outfile,'ptlevel');
fprintf('%s: Saved %s\n',mfilename,outfile);

metadata=ptlevel;

end

