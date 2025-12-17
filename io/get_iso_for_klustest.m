function [isods,quals,feats] = get_iso_for_klustest(dataformat,config,tets,clus,data_dirs)
% get_iso_for_klustest load spike features and calculate cluster quality
% Data loading function for klustest, given a data format type, and some
% settings, loads spike features calculated by kwiktint/klustakwik and
% calculates isolation distance and lratio.
%
% USAGE
%
% [isods,quals,feats] = get_iso_for_klustest(dataformat,config,tets,clus,data_dirs)
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
% 'clus' - Vector, Nx1 where N = number of spikes, gives the cluster assigned to
%       each spike, given by get_clu_for_klustest
%
% 'data_dirs' - Data directories, output from get_tets_for_klustest
%
% OUTPUT
%
% 'isods' - Cell array, isolation distance for each cluster on each tetrode,
%           cells correspond to electrodes
%
% 'quals' - Cell array, feature space histograms
%           (bin x-values; within cluster distances; between cluster distances)
%           For every cluster, cells correspond to electrodes
%
% 'feats' - Cell array, spike features for all spikes, cells correspond to
%           electrodes
%
% NOTES
% 1. Phy data format is experimental
%
% 2. Function is not really intended to be used without klustest
% 
% SEE ALSO kwiktint klustest mahal

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
        case {'kwikcut','kwiktint','phy','klustakwik'}
            isods = cell(1,size(tets,1));
            feats = cell(1,size(tets,1));
            quals = cell(1,size(tets,1));  

            for tt = 1:size(tets,1)
                if tt==1
                    disp(sprintf('\t\t...electrode: %d',tets(tt,1)))             
                else
                    disp(sprintf('\b, %d',tets(tt,1))) 
                end
                
                switch dataformat
                    case {'klustakwik','kwikcut','kwiktint'}
                        % get feature data                    
                        fetfile = [config.cname '.fet.' num2str(tets(tt,1))];
                        fid = fopen(fetfile,'r');
                        fets = sscanf(fgetl(fid),'%d');
    
                        if tets(tt,2)==1 % stereotrode
                            nch = 2;
                        elseif tets(tt,2)==0 % tetrode
                            nch = 4;
                        else % unknown electrode type
                            % get_tets_for_klustest listed an unknown electrode type here so we don't know how many channels it has
                            keyboard
                        end
    
                        if any(strcmp(dataformat,{'kwiktint'}))
                            nfets = (fets-1)/nch; % kwiktint fet files include spike times as the last column
                        else
                            nfets = fets/nch;
                        end
    
                        if round(nfets)~=nfets
                            % for some reason the number of features is not a whole number
                            keyboard
                        end
                        fdat = fscanf(fid,'%f',[fets Inf]);
                        fdata = fdat';
                        if any(strcmp(dataformat,{'kwiktint'}))
                            fdata = fdata(:,1:end-1); % kwiktint fet files include spike times as the last column
                        end                    
                        fetindx = max(fdata)~=min(fdata); % logical vector where 0 = dead channels
    
                        % get cluster data
                        clu = clus{1,tets(tt)}; 
                        fclose(fid);
    
                    case {'phy'}
                        fdata = [];
                        clu = [];
                        for ff = 1:length(data_dirs)
                            fetfile = fullfile(data_dirs{ff},[data_dirs{ff} '.mountainsort'],['output_T' num2str(tets(tt))],'phy_MS','pc_features.npy');
                            nfets = 5;  
                            nch = 4;
                            fnow = readNPY(fetfile); 
                            fnow = [fnow(:,:,1) fnow(:,:,2) fnow(:,:,3) fnow(:,:,4) zeros(size(fnow,1),1)];
                            fdata = [fdata; fnow];
    
                            cutname = fullfile(data_dirs{ff},[data_dirs{ff} '.mountainsort'],['output_T' num2str(tets(tt))],'phy_MS','spike_clusters.npy');
                            cnow = readNPY(cutname);  
                            clu = [clu; cnow];
                        end
                        fetindx = max(fdata)~=min(fdata); % logical vector where 0 = dead channels
                        
                    otherwise
                        keyboard
                        
                end
                clusters = unique(clu(clu>0));            
                isod = NaN(numel(clusters),1); % prepare arrays for data
                lrat = NaN(numel(clusters),1);
                dists = cell(max(clusters),1);
    
                %% run through each cluster and determine its properties
                for cc = 1:numel(clusters)
                    cnow = clusters(cc); % current cluster
                    % skip the noise cluster
                    if ~cnow
                        continue
                    end
                    cindx = clu == cnow; % index into .cut file for current cluster
    
                    % determine spikes inside vs outside cluster
                    spikes_in = sum(cindx); % number of cluster spikes
    
                    % Compute Mahalanobis distance between spikes inside (mdist_in) and outside (mdist_out) the cluster
                    mdist = mahal(fdata(:,fetindx),fdata(cindx,fetindx));
                    mdist_in = mdist(cindx);
                    mdist_out = mdist(~cindx);
    
                    % generate histograms of Mahalanobis distance for plottng later
                    zcut = 5;            
                    i = (mdist - mean(mdist(:))) ./ std(mdist(:)); % to remove spurious values which ruin the final plot            
                    b = logspace(log10( min(mdist(i<zcut)) ), log10( max(mdist(i<zcut)) ),128);
                    xi = movmean(b,2,'EndPoints','discard');
                    y1 = histcounts(mdist_in,b,'Normalization','count');
                    y2 = histcounts(mdist_out,b,'Normalization','count');
                    dists{cnow} = single([xi;y1;y2]);
    
                    % figure
                    % d = double(dists{cnow});
                    % y1 = d(2,:);%./sum(d(2,:));
                    % y2 = d(3,:);%./sum(d(3,:));
                    % b1 = bar(d(1,:),y1,1,'b','FaceAlpha',0.5,'EdgeColor','none'); hold on;
                    % b2 = bar(d(1,:),y2,1,'k','FaceAlpha',0.5,'EdgeColor','none');   
                    % ax = gca;
                    % ax.XScale = 'log';
    
                    % we need at least as many spikes as features, but not more than half the total number of spikes           
                    if spikes_in < (nfets*nch) || spikes_in > numel(clu)/2 
                        isod(cnow) = NaN;
                        lrat(cnow) = NaN;
                    else
                        % Determine the Mahalanobis of the Nth closest spike outside the cluster (where N is the number of spikes inside the cluster)
                        mdist_out_sorted = sort(mdist_out,'ascend');
                        isod(cnow) = mdist_out_sorted(spikes_in);    
        
                        % Calculate l-ratio
                        lquant = sum(1 - chi2cdf(mdist_out,nfets)); % the L quantity
                        lrat(cnow) = lquant / spikes_in; % Lratio is defined as L divided by the total number of spikes in the cluster                        
                    end
                end
                isods(1,tets(tt)) = {[isod(:) lrat(:)]}; % [isolation distance, lratio]
                feats(1,tets(tt)) = {int16(fdata)}; % feature space, we don't need a lot of precision for plotting
                quals(1,tets(tt)) = {dists}; % [bin x-values; within cluster distances; between cluster ditances]  
            end
    end




















































