function [clus] = get_clu_for_klustest(dataformat,config,tets,data_dirs)
% get_clu_for_klustest load clusters from .cut files or equivalent
% Data loading function for klustest, given a data format type, and some
% settings, finds .cut files (Tint format) or .cut files (Neuralynx format) and
% extracts cluster definitions.
%
% USAGE
%
% [clus] = get_clu_for_klustest(dataformat,config,tets,data_dirs)
%
% INPUT
%
% 'dataformat' - String: 'kwikcut','neuralynx','kwiktint','phy' (experimental)
%
% 'config' - Structure, configuration file from klustest, contains the
%           input/output filename to check for tetrodes etc
%
% 'tets' - Vector, electrodes we want to analyse
%
% 'data_dirs' - Data directories, output from get_tets_for_klustest
%
% OUTPUT
%
% 'clus' - Vector, Nx1 where N = number of spikes, gives the cluster assigned to
%       each spike.
%
% NOTES
% 1. Phy data format is experimental
%
% 2. Function is not really intended to be used without klustest
% 
% SEE ALSO kwiktint klustest getcut

% HISTORY
%
% version 1.0.0, Release 24/08/16 Initial release, decided to create specific loading functions for klustest
% version 1.0.1, Release 04/04/18 Commenting and formatting for klustest update
% version 1.0.2, Release 29/11/25 Updates for GitHub release
% version 1.0.3, Release 16/12/25 updated filenames for cross platform flexibility
% version 1.0.4, Release 16/12/25 improved comments, not ideal but detailed enough for now
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2025 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
switch dataformat
    case {'kwiktint'}
        clus = cell(1,size(tets,1));
       
        for tt = 1:size(tets,1)
            [clu,~] = getcut([config.cname '_' num2str(tets(tt)) '.cut']); % this reads the .cut file line by line so that it can be entered into the mtint structure  
            clus{tt} = uint8( clu(:) );  

            nclus = numel(unique(clus{tt}));
            disp(sprintf('\t\t...%d total clusters %s',nclus))
        end

    case {'kwikcut'}
        clus = cell(1,size(tets,1));
       
        for tt = 1:size(tets,1)
            electrode_type = tets(tt,2);
            if electrode_type==0 % tetrode
                cutname = [config.cname '_' num2str(tets(tt,1)) '.ntt']; % the cut file name                
            elseif electrode_type==1 % stereotrode
                cutname = [config.cname '_' num2str(tets(tt,1)) '.nst']; % the cut file name                                
            end
            clu = Nlx2MatSpike(cutname,[0 0 1 0 0],0,1,1); 
            clus{tt} = uint8( clu(:) );  

            nclus = numel(unique(clus{tt}));
            disp(sprintf('\t\t...%d total clusters %s',nclus))
        end
        
    case {'neuralynx'}
        clus = cell(1,size(tets,1));
        for tt = 1:size(tets,1)
            cutname = [config.cname '_' num2str(tets(tt,1)) '.cut']; % the cut file name
            clus{tt} = uint8( getcut(cutname) );  

            nclus = numel(unique(clus{tt}));
            disp(sprintf('\t\t...%d total clusters %s',nclus))
        end
        
    case {'phy'}
        clus = cell(1,max(tets));
        for tt = 1:length(tets)
            for ff = 1:length(data_dirs)
                cutname = fullfile(data_dirs{ff},[data_dirs{ff} '.mountainsort'],['output_T' num2str(tets(tt))],'phy_MS','spike_clusters.npy');
                clu = readNPY(cutname);

                % ignore noise clusters
                fname = fullfile(data_dirs{ff},[data_dirs{ff} '.mountainsort'],['output_T' num2str(tets(tt))],'phy_MS','cluster_group.tsv');
                fid = fopen(fname);
                A1 = textscan(fid,'%s\t%s',1); % column headers
                A2 = textscan(fid,'%d\t%s'); % {clusters} {groups}
                fclose(fid);     
                clu_good = A2{1,1}(ismember(A2{1,2},{'good'}));
                clu(~ismember(clu,clu_good)) = 0;
                clus{tets(tt)} = [clus{tets(tt)}; uint8( clu(:) )];  
            end
            clu_n = unique(clu);
            nclus = length(clu_n);
            disp(sprintf('\t\t...%d total clusters %s',nclus))
        end

end



































