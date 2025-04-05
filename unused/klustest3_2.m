function klustest3_2(cname,tetrodes,clusters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function utilises many minor functions to both analyse and plot cluster data generated 
%   in Tint. If requested, it will output a figure for each cluster, the cluster space of each
%   tetrode, the cross-correlations of ever cluster on a tetrode and a session data structure: sdata.mat.
%   It will also generate an mtint file (mtint.mat) containing all the tetrode and cluster info.
%   KlusterAnalysisTINT(tetrodes,clusters)
%
%   3D tracking data can be edited to remove outliers using editLATTICEpos
%
%%%%%%%% Inputs
%   tetrodes = (default = 1:16) the tetrodes to run on in a vector format (i.e. [1 2 3 4 5 6] would run on tetrodes 1 to 6)
%   clusters = (default = 0) the clusters to run on, set this to 0 to run on all clusters
%
%%%%%%%% Comments
%   05/08/16 created from an older version
%   08/08/16 modified readDACQDATA to output waveforms and added post processing of this info
%   08/08/16 readDACQDATA working, added postprocessing of mtint
%   09/08/16 main plots working
%   10/08/16 options added for different styles of plot
%   11/08/16 added theta autocorrelation sine estimate (O'Mara 2015)
%   12/08/16 added refractory period violations (Navratilova and McNaughton (2016)
%   13/08/16 added saving data, figures
%   14/08/16 added cell type identification
%   15/08/16 cluster quality assessment added in clusterQUALITY
%   16/08/16 added cluster space figures
%   17/08/16 fixed phase plot and phase map
%   18/08/16 modified clusterQUALITY to deal better with missing channels
%   18/08/16 added the option to ignore position data, fixed bug with waveform plot
%   19/08/16 make it so none of the position figures are made if there is no position data, this should be faster
%   22/08/16 fixed issues with cluster space plot, concentrate on first feature, added legend
%   23/08/16 added cluster space subplot to cell figure, for Ele, this uses the first clustering feature, the 1st and 2nd highest amplitude channel
%   25/06/16 added tetrode input to readalldacqdata, this prevents it trying to open unnecessary files, fixed bug in the detection of dead/lfp channels
%   25/08/16 fixed a minor bug in histograms when there is only one spike
%   22/11/16 adapted KlusterAnalysisTINT for 3D mazes, renamed klutest3
%   16/02/17 removed unnecessary settings
%   24/03/17 klustest3_1 created
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manual settings
close all
fclose('all');

% initial variables
rname = '851';

% overrides
mtint_override = 0; % set this to 1 to force the production of a new mtint file (do this if you reclustered in tint for instance)
frame_override = 0; % set to 1 to override any saved lattice/frame 
expand_override = 0; % set to 1 to recompute dwellmap/expansion

% partitioning
part_config{1} = {'square1',2,1,[],2}; % part name, method, intervals, times, dimensions
part_config{2} = {'lattice',2,2,[],3};
part_config{3} = {'square2',2,3,[],2};
% part_config{1} = {'lattice',2,2,[],3};

% 3d settings
config3 = struct;
config3.pratio = 1000;
config3.leds = 1;
config3.srate = 50;
config3.lat_size = [970 970 970]; % lattice size in mm
config3.of_size = [1200 1200 200]; % open field size in mm
config3.vsize = [25 25 25]; % voxel size in mm [x,y,z] or [v]
config3.sigma = [50 50]; % standard deviation for position and spike data gaussian [p s] or [x]
config3.mindist = 100; % minimum distance (mm) from position data for a voxel to be filled
config3.mindwell = 1; % minimum time (s) for considering a voxel visited
config3.psize = [4 4 4]; % number of padding voxels [x,y,z] or [v]

% place field settings
frcut = 0.2; % percentage cutoff for ratemaps, before detecting place fields
vlcut = 9^3; % the minimum volume (bins) of a place field   

% figure settings
fig_vis = 'off';
save_fig = 1; % set to 1 to save .fig files for 3D viewing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% If no tetrodes are specified assume all of them (these will be checked later)
if ~exist('tetrodes','var') || isempty(tetrodes)
    tetrodes = 1:16;
end % if ~exist('tetrodes','var') || ismepty(tetrodes)

%% If no clusters are specified assume all of them
if ~exist('clusters','var') || isempty(clusters)
    clusters = 0;
end % if ~exist('clusters','var') || ismepty(clusters)

%% If no output name is specified, assume kwiktint output files with name 'kwiktint'
if ~exist('cname','var') || isempty(cname)
    cname = 'kwiktint';
end % if ~exist('clusters','var') || ismepty(clusters)

%% Start analysis
tic;
disp('----------------------------------------------------------------------------');
disp(sprintf('Running klustest...'))
sdata = struct; % create an empty structure - together with the mtint file this will hold all of the session data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get session data
disp(sprintf('Identifying sessions...'))
[snames,~,nsess] = getSNAMES;
sdata.session_names = snames; % add data to structure
sdata.combined_name = cname; % add data to structure
sdata.num_sessions = nsess;

%% Create directories								                
disp('Preparing directories...')
[~,~,~] = mkdir([pwd '\klustest3\figures']);
disp(sprintf('\t...done'))

%% Check tetrodes
disp(sprintf('Assessing data...'))
[tetrodes,mvalue] = getTRODES(snames,tetrodes);
if isempty(mvalue)
    disp(sprintf('\t...tetrodes: %s accounted for',mat2str(tetrodes)))
else
    disp(sprintf('\t...WARNING: tetrodes: %s are incomplete and will be skipped',mat2str(mvalue)))
end % if isempty(mvalue)
disp(sprintf('\t...done'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read all dacq data or load it
disp(sprintf('Fetching DACQ data...'));
if any(strcmp(evalin('base','who'),'mtint')) && ~mtint_override
    disp(sprintf('\t...using mtint held in memory'));
    mtint = evalin('base','mtint');
elseif ~exist(['klustest3\' cname '_mtint.mat'],'file') || mtint_override
    disp(sprintf('\t...running getDACQDATA'));
    mtint = getDACQDATA(cname,snames);
    disp(sprintf('\t...post-processing mtint'));
    mtint = postprocess_DACQ_data(mtint);
    info = whos('mtint');
    siz = info.bytes / 1000000;
    disp(sprintf('\t...saving %.1fMb mtint',siz))
    save(['klustest3\' cname '_mtint.mat'],'cname','mtint','-v7.3');
else
    disp(sprintf('\t...loading saved mtint'));
    load(['klustest3\' cname '_mtint.mat'],'-mat');
end % if ~exist(['klustest3\' cname '_mtint.mat'],'file') || mtint_override
disp(sprintf('\t...done'));
sdata.date = mtint.header_date;

assignin('base','mtint',mtint); % leave mtint in base workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get some variables that will be important later
disp(sprintf('Extracting initial data...'))
disp(sprintf('\t...recording sessions: %d',numel(snames)))

%% Get the position data for the sessions
% manually define some settings
samp_rate_hz = config3.srate; % synchronised dacqTrack data should be 50Hz
pos_tb = samp_rate_hz / 1000;    
pixel_ratio = config3.pratio;
leds = config3.leds; % typically we only use 1 LED

pox = []; poy = []; poz = []; % initialise variables
for pp = 1:nsess % for each session
    snow = snames{pp}; % get its filename
    fcheck = [pwd '\triangulation\' snow '_merge_data.mat']; % load its 3d reconstructed path
    if exist(fcheck,'file') % check it exists
        load(fcheck);
    else
        error('Predetermined 3D trajectory not found: %s\n ...exiting',fcheck)
    end % if exist(fcheck,'file')

    posn = datout.weight_mean_merge_smoothed; % extract the smoothed, weighted mean trajectory
    scount_now = mtint.pos.trial_samples(pp); % only take the first n samples, where n = the number of samples in the dacqUSB data (often there is one extra data point in the synchronised data)
    dcount = abs(length(posn(:,1))-scount_now);
    if length(posn(:,1)) < scount_now
        disp(sprintf('\t...3D session shorter by %d points (%.2fs), NaN padding data',dcount,(dcount * (1/samp_rate_hz))));
        posn2 = NaN(scount_now,3);
        posn2(1:size(posn,1),1:size(posn,2)) = posn;
        posn = posn2;
    elseif length(posn(:,1)) > scount_now
        disp(sprintf('\t...3D session longer by %d points (%.2fs), trimming data',dcount,(dcount * (1/samp_rate_hz))));   
        posn = posn(1:scount_now,:);        
    end % if length(posn(:,1)) < scount_now

    pox = [pox; posn(:,1)]; % extract the x data and concatenate it
    poy = [poy; posn(:,2)]; % extract the y data and concatenate it
    poz = [poz; posn(:,3)]; % extract the z data and concatenate it
end % for ss = 1:nsess
pox = pox(:);
poy = poy(:);
poz = poz(:);

poss1 = numel(mtint.pos.xy_pixels(:,1));
poss2 = numel(pox);
disp(sprintf('\t...total 3D tracking positions read: %d',poss2));
if poss1 ~= poss2
    error('ERROR: different number of data points in 2D (%.f) and 3D (%.f) sessions... exiting',poss1,poss2);
end % if poss1 ~= poss2
pot = double(mtint.pos.ts); % the 3D pos times should be the same as the regular data
sdata.pox = single(pox); % add data to structure
sdata.poy = single(poy); % add data to structure
sdata.poz = single(poz); % add data to structure
sdata.pot = single(pot); % add data to structure
head_direct = zeros(size(pox));

% get the total length of the session
duration = key_value('duration',mtint.pos.header,'num');
disp(sprintf('\t...total session time: %ds',duration))
sdata.session_duration = duration; % add data to structure

% display results
disp(sprintf('\t...positions read: %d',poss1));
disp(sprintf('\t...tracking LEDs: %d',leds));
disp(sprintf('\t...pixel ratio: %dppm',pixel_ratio));
disp(sprintf('\t...sample rate: %dHz (%.2fs)',samp_rate_hz,pos_tb));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve and analyse data
disp(sprintf('Analysing lfp data...'))
lfp = double(mtint.lfp(1).lfp(:,1));
Fs = mtint.lfp(1).Fs(1,1);
lfptime = (0:length(lfp)-1)'/Fs; % make a vector for time
if round(max(lfptime)) ~= round(duration)
    disp(sprintf('WARNING: lfp time %.2f does not match session duration %.2f...',max(lfptime),duration));
end % if max(lfptime) ~= duration

% filter lfp to get theta
cutOffFreq = [4 12]; % lower and upper cutoff frequencies in Hz
[b,a] = butter(4,cutOffFreq/(Fs/2)); % Generate 4th order butterworth filter coefficients
lfpfilt = filtfilt(b,a,lfp); % Apply filter to data using zero-phase filtering
disp(sprintf('\t...done'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% partitioning settings
disp(sprintf('Preparing partitions...'))
part_config = partSESS(part_config,mtint);
nparts = numel(part_config);
disp(sprintf('\t...done'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve and analyse data
disp(sprintf('Analysing tetrode data...'))
%% For every available tetrode
for e = 1:length(tetrodes) 
    tet = tetrodes(e); % tet = the current tetrode
    disp(sprintf('\tLoading electrode %d...',tet));
    
    pos_samps = ceil(mtint.tetrode(tet).ts * 50);    
    spik_count = numel(mtint.tetrode(tet).nspike_cut); % retrieve the number of spikes recorded on this tetrode

    % get a vector of clusters we want to analyse
    if clusters == 0
        clus = unique(mtint.tetrode(tet).cut); % find the clusters logged for this tetrode
        clus(clus == 0) = []; % find and delete cluster 0 
    else
        clus = clusters;
    end % if clusters == 0
    
    % check to see if there are any clusters
    clus_count = numel(clus);   
    if ~clus_count % if there are no clusters
        continue % skip analysis
    end % if ~clus_count

    disp(sprintf('\t\t...%d spikes detected',spik_count));
    disp(sprintf('\t\t\t...%d data clusters detected',clus_count));    
    disp(sprintf('\t\t\t\t...starting analysis'));
    disp(sprintf('\t\t\t\t...progress: 0%%'));

    % get the spike times for this tetrode
    spiketime = mtint.tetrode(tet).ts;

    % calculate the quality of clusters in this session
    [isods,lratios,~,~,~,~] = clusterQUALITY(cname,tet,clus);        
        
    % get the channel waveforms for this tetrode
    waves = cell(4);
    for ggnow = 1:length(sdata.session_names)
        fnamen = sdata.session_names{ggnow};
        [~,c1,c2,c3,c4] = getspikes([fnamen,'.',num2str(tet)]);   
        waves{1} = [waves{1}; c1];
        waves{2} = [waves{2}; c2];
        waves{3} = [waves{3}; c3];
        waves{4} = [waves{4}; c4];  
    end % for ggnow = 1:length(sdata.session_names)      
    
%% For every detected cluster           
    for cc = 1:length(clus)
        clu = clus(cc); % clu = the current cluster
        uci = ['u' rname sdata.date num2str(tet) num2str(clu)]; % unique cell identifier - a string with [rat name, session date, electrode, cluster];
            
        isod = isods(cc); % retrieve cluster quality info
        lratio = lratios(cc);
        sdata.(uci).cluster_isolation = isod; % add data to structure
        sdata.(uci).cluster_lratio = lratio; % add data to structure
           
        % get some vectors that we can use to sort data
        clu_identity = mtint.tetrode(tet).cut; % clu_assign is a vector of numbers, one for each spike, each number corresponds to a cluster
        clu_indx = (clu_identity == clu); % logical vector, one number per spike, 1 if in cluster, 0 if not
        n_spikes = length(find(clu_identity == clu));
        pos_identity = pos_samps; % pos_identity is a vector of numbers, one for each spike, each number corresponds to a position data point
            
        % sort out cluster spike data
        cspt = spiketime(clu_indx); % the time point for every spike in this cluster        
        spkindx = knnsearch(pot,cspt); % get an index of position data for each spike (i.e. which pox,poy,poz value most closely corresponds to each spike)    
        cspx = pox(spkindx);
        cspy = poy(spkindx);                       
        cspz = poz(spkindx);    
        cspx = cspx(:);
        cspy = cspy(:);
        cspz = cspz(:);

        %sdata.(uci).spx = single(cspx); % add data to structure
        %sdata.(uci).spy = single(cspy); % add data to structure
        %sdata.(uci).spz = single(cspz); % add data to structure            
        %sdata.(uci).spt = single(cspt); % add data to structure                      
        chd = double(head_direct(pos_identity(clu_identity == clu))); % get an index of which spikes belong to this cluster, then get an index of whcih positions these spikes correspond to, then use this to get the correct head directions       

%% For each part       
        for pp = 1:nparts % for every partition
            part_now = part_config{pp}{1}; % the name of the current part
            part_times = part_config{pp}{4}; % the time pairs (intervals) corresponding to this part
                
            ispk = []; % holder for interval spikes
            ipos = []; % holder for interval position data           
            part_duration = 0; 
            
%% For each interval   
            sindax = [];
            for ii = 1:length(part_times(:,1)) % for every pair of time values (interval) associated with this part  
                i_time = part_times(ii,:); % get the start and end time
                sindx = find(cspt > i_time(1) & cspt < i_time(2)); % find the spikes falling into this interval
                pindx = find(pot > i_time(1) & pot < i_time(2)); % find the position data falling into this interval

                % concatenate all the data
                ispk = [ispk; cspx(sindx) cspy(sindx) cspz(sindx) cspt(sindx) chd(sindx) ones(size(cspx(sindx))).*ii];
                ipos = [ipos; pox(pindx) poy(pindx) poz(pindx) pot(pindx) ones(size(pox(pindx))).*ii];                
                part_duration = part_duration + (i_time(2)-i_time(1));
                sindax = [sindax; sindx(:)]; % vector for pulling out spikes for this part later    
            end % for tt = 1:length(goal_times(:,1))   

            % assign data to arrays
            ppox = ipos(:,1);
            ppoy = ipos(:,2);
            ppoz = ipos(:,3);                
            ppot = ipos(:,4);
            pspx = ispk(:,1);
            pspy = ispk(:,2);
            pspz = ispk(:,3);                
            pspt = ispk(:,4);    
            sdata.(part_now).pox = ppox;
            sdata.(part_now).poy = ppoy;
            sdata.(part_now).poz = ppoz;
            sdata.(part_now).pot = ppot;            
            sdata.(uci).(part_now).spx = pspx;
            sdata.(uci).(part_now).spy = pspy;
            sdata.(uci).(part_now).spz = pspz;

            % work out the firing rate of the cluster in this part
            pfrate = numel(pspx) / part_duration;
            sdata.(uci).(part_now).frate = pfrate;
            sdata.(part_now).duration = part_duration;            
            
            % fit a frame to this partition data
            lname = ['klustest3\' part_now '_frame.mat'];
            if exist(lname,'file') && ~frame_override
                load(lname);
            else
                if part_config{pp}{5} == 2 % if this part is 2 dimensional
                    [lx,ly,lz] = makeLATTICE_v2(ppox,ppoy,ppoz,config3.pratio,config3.of_size,1); % we want a basic cuboid wire frame
                else
                    [lx,ly,lz] = makeLATTICE_v2(ppox,ppoy,ppoz,config3.pratio,config3.lat_size); % we want a 3D cube lattice
                end % if part_config{pp}{5} == 2
                save(lname,'lx','ly','lz') 
            end % if exist('3D matfiles\lattice.mat')
            sdata.(part_now).lx = lx;
            sdata.(part_now).ly = ly;
            sdata.(part_now).lz = lz;
                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Waveform data
            if ~numel(pspx) % if there are no spikes  
                continue
            end % if ~numel(gspx) % if there are no spikes
            
            % get the waveforms for this cluster
            waves2{1} = waves{1}(clu_indx,:);
            waves2{2} = waves{2}(clu_indx,:);
            waves2{3} = waves{3}(clu_indx,:);
            waves2{4} = waves{4}(clu_indx,:);                    
            waves2{1} = waves2{1}(sindax,:);
            waves2{2} = waves2{2}(sindax,:);
            waves2{3} = waves2{3}(sindax,:);
            waves2{4} = waves2{4}(sindax,:);
            for w = 1:4 % for every recording channel
                wav = double(waves2{w});
                ch = squeeze(mean(wav,1));
                chs = squeeze(std(wav,1));
                [maxval,maxt] = max(ch);
                [postminval,postmint] = min(ch(maxt:end));
                postmint = postmint + maxt - 1;
                width = postmint - maxt;
                width = width * (1000/50); 
                
                sdata.(uci).(part_now).waveform_mean{w} = ch;
                sdata.(uci).(part_now).waveform_stdv{w} = chs;
                sdata.(uci).(part_now).waveform_max(w) = maxval;
                sdata.(uci).(part_now).waveform_min(w) = postminval;
                sdata.(uci).(part_now).waveform_width(w) = width;
            end % for w = 1:4 % for every recording channel             
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate 3D ratemap and dwellmap               
            % retrieve a precomputed dwellmap if one exists
            dname = ['klustest3\' part_now '_dwellmap.mat'];
            if exist(dname,'file') && ~expand_override % if a precomputed dwellmap exists        
                load(dname,'dwellmap','config3'); % load in data
                [ratemap,dwellmap] = generate3Dmap_v2(dwellmap,[ppox ppoy ppoz],[pspx pspy pspz],config3.vsize,config3.sigma,config3.mindist,config3.mindwell,config3.srate,config3.pratio,config3.psize);         
            else
                [ratemap,dwellmap] = generate3Dmap_v2([],[ppox ppoy ppoz],[pspx pspy pspz],config3.vsize,config3.sigma,config3.mindist,config3.mindwell,config3.srate,config3.pratio,config3.psize);            
                save(dname,'dwellmap','config3'); % save data for later use
            end % if exist(dname,'file') 
            sdata.(uci).(part_now).ratemapA = single(ratemap);
            sdata.(part_now).dwellmapA = single(dwellmap); 
            
            % convert lattice vertices to match the ratemap   
            lx2 = lx ./ config3.vsize(1);     
            ly2 = ly ./ config3.vsize(2);                 
            lz2 = lz ./ config3.vsize(3);  
            lx2 = lx2 - mean([min(lx2) max(lx2)]);
            ly2 = ly2 - mean([min(ly2) max(ly2)]);
            lz2 = lz2 - mean([min(lz2) max(lz2)]);
            lx2 = lx2 + (size(ratemap,1)/2);
            ly2 = ly2 + (size(ratemap,2)/2);
            lz2 = lz2 + (size(ratemap,3)/2);    
            sdata.(part_now).lxrm = lx2;
            sdata.(part_now).lyrm = ly2;            
            sdata.(part_now).lzrm = lz2;
            
            % threshold ratemap
            %thresh = nanmean(nanmean(nanmean(ratemap)));
            thresh = frcut;
            ratemap1 = ratemap;
            dwellmap1 = dwellmap;
            ratemap(ratemap < thresh*nanmax(ratemap(:))) = NaN;                         
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%% Initial ratemap analyses       
            % spatial information content and sparsity
            spati = Skaggs3D(ratemap1,dwellmap1);
            spars = Spars3D(ratemap1,dwellmap1);     
            setpy = entropy(ratemap1);
            sdata.(uci).(part_now).sparsity = spars;
            sdata.(uci).(part_now).spatinfo = spati;
            sdata.(uci).(part_now).sentropy = setpy;

            % find place fields
            [fields,fdata] = fields3D(ratemap,frcut,vlcut,config3.vsize(1));                      
            sdata.(uci).(part_now).nfields = fields;
            sdata.(uci).(part_now).field_data = fdata;

            % skaggs compression matrices
            [skag1,skag2,skag3,mmap1,mmap2,mmap3] = SkaggsCompress3D(ratemap1,dwellmap1);
            sdata.(uci).(part_now).spatinfo_dims = [skag1,skag2,skag3];
            sdata.(uci).(part_now).spatinfo_maps = {mmap1,mmap2,mmap3};

            % compression matrix grid scores
            auto1 = xPearson(mmap1);
            auto2 = xPearson(mmap2);
            auto3 = xPearson(mmap3);
            [gscor1,gspac1,fs1,go1,ge1] = GridAnalysis(auto1,config3.vsize(1));      
            [gscor2,gspac2,fs2,go2,ge2] = GridAnalysis(auto2,config3.vsize(1));      
            [gscor3,gspac3,fs3,go3,ge3] = GridAnalysis(auto3,config3.vsize(1));  
            sdata.(uci).(part_now).gscore_comp = [gscor1,gscor2,gscor3]; 
            sdata.(uci).(part_now).gspacing_comp = [gspac1,gspac2,gspac3]; 
            sdata.(uci).(part_now).gfsize_comp = [fs1,fs2,fs3]; 
            sdata.(uci).(part_now).gorientation_comp = [go1,go2,go3]; 
            sdata.(uci).(part_now).gellipticity_comp = [ge1,ge2,ge3]; 
            sdata.(uci).(part_now).gridscore_maps = {auto1,auto2,auto3};

            % map slices
            xs = round(unique(lx2));
            xs = xs(~isnan(xs));
            xs(xs==0) = 1;
            ys = round(unique(ly2));
            ys = ys(~isnan(ys));
            ys(ys==0) = 1;                    
            zs = round(unique(lz2));
            zs = zs(~isnan(zs));
            zs(zs==0) = 1;                    

            [gx,gy,gz,gxm,gym,gzm,gxa,gya,gza] = GridAutoSlice3D(ratemap1,config3.vsize(1),xs,ys,zs);
            [sx,sy,sz,sxm,sym,szm,sxd,syd,szd] = SkaggsAutoSlice3D(ratemap1,dwellmap1,xs,ys,zs);                    
            sdata.(uci).(part_now).gscore_slice = {gx,gy,gz}; 
            sdata.(uci).(part_now).spatinfo_slice = {sx,sy,sz}; 
            sdata.(uci).(part_now).slices = {sxm,sym,szm}; 
            sdata.(uci).(part_now).auto_slices = {gxa,gya,gza}; 

            % spike distributions in lattice maze
            pos = [ppox ppoy ppoz ppot];                    
            spk = [pspx pspy pspz pspt];                  
            sdists = spikeDISTRIBUTION(spk,pos,{lx ly lz},0);
            sdata.(uci).(part_now).spike_dists = sdists; 
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures
            % Figure 1, ratemap, spike plot etc
            kfig_overall(sdata,uci,part_now,fig_vis,save_fig); % overall figure with ratemap etc

            % Figure 2 & 3, ratemap slices etc
            kfig_slices(sdata,uci,part_now,1,fig_vis,save_fig); % normal slices
            kfig_slices(sdata,uci,part_now,2,fig_vis,save_fig); % autocorrelation slices

        end % for gg = 1:length(part_indx) % for each specified goal
        disp(sprintf('\b %.f%%',cc/length(clus)*100));
    end % for c = 1:length(clusters) % for every detected cluster
end % for e = 1:length(tetrodes) % for every available tetrode         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Position data analysis
latticePATH(sdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finish up
% Save the session data structure
info = whos('sdata');
siz = info.bytes / 1000000;
disp(sprintf('\t...saving %.1fMb sdata',siz))
save(['klustest3\' cname '_sdata.mat'],'sdata','-v7.3'); % save session data

% final messages
toc1 = toc/60;
disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd '\Figures',' &'');">','Figures folder','</a>'])
disp(['Open ','<a href = "matlab: [s,r] = system(tint_location);">','TINT','</a>'])
disp(sprintf('klustest3 has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
disp(sprintf('Done.'));
disp('-------------------------------------------------------------------------------------');
































