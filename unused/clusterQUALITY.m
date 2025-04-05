function [isods,lrats,nds,cds,fdata,nfets,clus,fdata_cut] = clusterQUALITY(fname,tet,clu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial variables
if ~exist('clu','var') || isempty(clu)
    clu = -1;
end % if ~exist('clu','var') || isempty(clu)

isod = NaN;
lrat = NaN;
cd = [];
nd = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% import feature data
% in these each channel's data is concatenated into adjacent columns, with each column representing a feature (of which there can be any number)
% i.e. in an example where we have selected three PCs, then the .fet file would consist of:
% 13 (the first line states the total number of features which in this case is 3 PCs * 4 channels + 1 time
% PC1_1_1 PC2_1_1 PC3_1_1 PC1_2_1 PC2_2_1 PC3_2_1 PC1_3_1 PC2_3_1 PC3_3_1 PC1_4_1 PC2_4_1 PC3_4_1 TIMESTAMP_1 (so the channels are arranged ch1 features, ch2 features etc then time)

% get feature data
fetfile = [pwd '\' fname '.fet.' num2str(tet)];
fid = fopen(fetfile,'r');
fets = sscanf(fgetl(fid),'%d');
nfets = (fets-1)/4;
fdat = fscanf(fid,'%f',[fets Inf]);
fdata = fdat';

% get cluster data
cutfile = [pwd '\' fname '_' num2str(tet) '.cut'];
clufile = [pwd '\' fname '.clu.' num2str(tet)];
if exist(cutfile,'file')
    [clus,~] = getcut(cutfile); % this reads the .cut file line by line
elseif exist(clufile,'file')
    fid = fopen(clufile,'r');
    clus = fscanf(fid,'%d');
else
    clus = zeros(length(fdata(:,1)),1);
end % if exist([pwd '\' fname '_' num2str(tet) '.cut'],'file')

%% remove feature data from missing channels
% open Tint output file
fID = fopen([pwd '\kwiktint\TINT_output_' fname '_' num2str(tet) '.txt'],'r');
tdata = textscan(fID,'%s','Delimiter','\r');
tdata = tdata{:};
useindx = strfind(tdata,'UseFeatures'); % find lines containing this text

index = false(1,numel(useindx)); % create a logical array
for k = 1:numel(useindx) % for every line
  index(k) = isempty(useindx{k}); % check if it is empty
end % for k = 1:numel(sindx)
fetline = find(~index,1,'last'); % the last line listing features must have been the one used for the current analysis
fetline = tdata{fetline}; % extract line of text
secindx = strfind(fetline,'UseFeatures'); % find where it says what features were used
spcindx = strfind(fetline,' '); % also find all spaces in this line

sindx1 = find(spcindx > secindx,1,'first'); % the feature string starts with the first space after it says UseFeatures
srt_fet = spcindx(sindx1)+1; % this is the start of the string 
end_fet = spcindx(find(spcindx > srt_fet,1,'first'))-1; % the feature string ends with the next space
feature_string = fetline(srt_fet:end_fet);

fets_s = zeros(1,length(feature_string)); 
fets_s(strfind(feature_string,'1')) = 1; % convert the string to a logical array
fdata_cut = fdata;
fdata_cut(:,~fets_s) = []; % remove the dead  from our feature matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare data
if clu == -1;
    clusters = unique(clus);
    clusters(clusters == 0) = [];
else
    clusters = clu;
end % if clu == -1;

%% run through each cluster and determine its properties
nclus = numel(clusters);
isods = NaN(nclus);
lrats = NaN(nclus);
nds = cell(nclus);
cds = cell(nclus);
for cc = 1:nclus
    cnow = clusters(cc);
    cindx = find(clus == cnow);
    nindx = find(clus ~= cnow);
    clustfets = double(fdata_cut(cindx,:));
    noisefets = double(fdata_cut(nindx,:));

    a_spikes = numel(clus); % the number of spikes on this tetrode
    n_spikes = numel(cindx); % the number of spikes in this cluster

    siz1 = size(noisefets);
    siz2 = size(clustfets);

    if siz1(2) >= siz1(1) % if there are more population spikes than cluster features we can't calculate quality
        isod = -1;
        lrat = -1;
        return
    end % if s1(2) >= s1(1)

    if siz2(2) >= siz2(1) % if there are more cluster features than cluster spikes we can't calculate quality
        isod = -1;
        lrat = -1;
        return
    end % if s1(2) >= s1(1)

    %% calculate mahalanobis distance
    nd = mahal(noisefets,clustfets);
    nd = sort(nd,'ascend');
    cd = mahal(clustfets,clustfets);
    cd = sort(cd,'ascend');
    nds{cc} = nd;
    cds{cc} = cd;

    if n_spikes <= (a_spikes - n_spikes)
        isod = nd(n_spikes);
    else
        isod = -1;
    end % if n_spikes <= (a_spikes - n_spikes)
    isods(cc) = isod;
    
    lrat = sum(1-chi2cdf(nd,length(clustfets(1,:)))) / n_spikes;
    lrats(cc) = lrat;
end % for cc = 1:nclus
    













