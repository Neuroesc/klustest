function fetdata = clusterQUALITY_v2(fetfile,cutfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Schmitzer-Torbert (2005) http://dx.doi.org/10.1016/j.neuroscience.2004.09.066
%     Isolation Distance:
%         Isolation Distance is therefore the radius of the smallest ellipsoid from the cluster center containing all of the cluster spikes 
%         and an equal number of noise spikes. As such, Isolation Distance estimates how distant the cluster spikes are from the other spikes 
%         recorded on the same electrode. Isolation Distance is not defined for cases in which the number of cluster spikes is greater than 
%         the number of noise spikes.
%     Lratio:
%         Noise spikes which are close to the center of cluster C will contribute strongly to this sum, while noise spikes far from the center 
%         of cluster C will contribute little. A low value of L indicates that the cluster has a good “moat” and is well separated from other 
%         spikes recorded on the same tetrode. In contrast, a high value of L indicates that the cluster is not well separated, and is likely 
%         to both include spikes which are not part of the cluster and exclude spikes that are part of the cluster.
%
%   This function uses Tint files (.fet, .fmask and .cut) to calculate cluster quality measures
%   isolation distance (isod) and lratio (lrat). It also provides the mahalanobis distance
%   values computed during the process which can be plotted to visualise the effect.
%   For this calculation, noise is taken as all spikes on a tetrode that are not part of the
%   cluster. If there are too few of these then quality cannot be measured and both measures will
%   equal -1.
%   [isods,lrats,nds,cds,fdata,nfets,clus,fdata_cut] = clusterQUALITY(fname,tet,clu)
%
%%%%%%%% Inputs
%   fname = the filename of the Tint files being processed
%   tet = the tetrode to analyse
%   clu = the cluster to analyse
%   fdata = the extracted feature data, masked, each group of four columns represent a feature, with each column representing a channel, the final column is time
%   nfets = the number of features detected
%   clus = a vector where each number corresponds to a cluster assignment
%   fetNum = the feature string for this tetrode
%   fdata_cut = fdata minus any missing channels
%
%%%%%%%% Outputs
%   isods = isolation distances, this measures the distance a cluster is separated from other clusters (higher values are better)
%   lrats = lratios, this measures the amount of contamination inside a cluster (lower values are better)
%   nds = the mahalanobis distances between all spikes and the noises
%   cds = the mahalobis distances between all spikes in the clusters
%
%%%%%%%% Comments
%   09/08/16 created 
%   25/11/16 importdata removed and replaced with freadf etc for speed
%   25/11/16 changed function so its able to process multiple clusters, again for speed
%   04/04/17 created from clusterQUALITY, aim to fix error in l-ratio calculation and remove reliance on kwiktint outputs
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% import feature data
% in these each channel's data is concatenated into adjacent columns, with each column representing a feature (of which there can be any number)
% i.e. in an example where we have selected three PCs, then the .fet file would consist of:
% 13 (the first line states the total number of features which in this case is 3 PCs * 4 channels + 1 time
% PC1_1_1 PC2_1_1 PC3_1_1 PC1_2_1 PC2_2_1 PC3_2_1 PC1_3_1 PC2_3_1 PC3_3_1 PC1_4_1 PC2_4_1 PC3_4_1 TIMESTAMP_1 (so the channels are arranged ch1 features, ch2 features etc then time)

% fetfile = 'kwiktint.fet.5';
% cutfile = 'kwiktint_5.cut';

% get feature data
fid = fopen(fetfile,'r');
fets = sscanf(fgetl(fid),'%d');
nfets = (fets-1)/4;
fdat = fscanf(fid,'%f',[fets Inf]);
fdata = fdat';
fetindx = max(fdata)~=min(fdata); % logical vector where 0 = dead channels

% get cluster data
clus = getcut(cutfile); % this reads the .cut file line by line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare data
clusters = unique(clus);
nclus = numel(clusters); % count the number of clusters
isods = NaN(max(clusters),1); % prepare arrays for data
lrats = NaN(max(clusters),1);
noise_dists = cell(max(clusters),1);
clust_dists = cell(max(clusters),1);

%% run through each cluster and determine its properties
for cc = 1:nclus
    cnow = clusters(cc); % current cluster
    % skip the noise cluster
    if ~cnow
        continue
    end
    cindx = clus == cnow; % index into .cut file for current cluster
    
    % determine spikes inside vs outside cluster
    tspikes = numel(clus); % total number of spikes
    cspikes = sum(cindx); % number of cluster spikes
    nspikes = tspikes - cspikes; % number of noise spikes (i.e. spikes not in the current cluster)
    
    % we need at least as many spikes as features, but not more than half the total number of spikes
    if cspikes < nfets || cspikes > tspikes/2 
        isods(cnow) = -1;
        lrats(cnow) = -1;
        continue
    end
    
    % Compute Mahalanobis distance between spikes inside (mdist_in) and outside (mdist_out) the cluster
    mdist = mahal(fdata(:,fetindx),fdata(cindx,fetindx));
    mdist_in = mdist(cindx);
    mdist_out = mdist(~cindx);
    noise_dists{cnow} = uint16(mdist_out);
    clust_dists{cnow} = uint16(mdist_in);

    % Determine the Mahalanobis of the Nth closest spike outside the cluster (where N is the number of spikes inside the cluster)
    mdist_out_sorted = sort(mdist_out,'ascend');
    isod = mdist_out_sorted(cspikes);    
    isods(cnow) = isod;
    
    % Calculate l-ratio
    lquant = sum(1 - chi2cdf(mdist_out,nfets)); % the L quantity
    lrat = lquant / cspikes; % Lratio is defined as L divided by the total number of spikes in the cluster
    lrats(cnow) = lrat;
    
end

%% accumulate data
fetdata = struct;
fetdata.fetdata = int16(fdata);
fetdata.nclusters = nclus;
fetdata.nspikes = length(fdata(:,1));
fetdata.nfeatures = nfets;
fetdata.fetindx = fetindx;
fetdata.isolation_distances = isods;
fetdata.lratios = lrats;
fetdata.noise_dists = noise_dists;
fetdata.clust_dists = clust_dists;
% fetdata.cluster_labels = clus;













