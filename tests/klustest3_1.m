function klustest3_1(tetrodes,clusters)
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

% overrides
mtint_override = 0; % set this to 1 to force the production of a new mtint file (do this if you reclustered in tint for instance)
elook_override = 0; % set to 1 to override an existing expansion lookup table (i.e. if you change any 3D ratemap settings this will have to be overwritten)

% partitioning
sess_data.names = {'square_a','lattice','square_b'}; % the names of the sessions we ultimately want to end up with
sess_data.assignment = [1,2,3]; % the 
sess_data.dimensions = [2,3,2];
sess_data.ignore = [0,0,0];

% 3d settings
ratio_3d = 1000; % pixels per metre, put 1000 for mm accurate 3D reconstruction
bsize_3d = 5; % starting binsize in cm3, can also be a 3 element vector i.e. [10 10 20] for anisotropic bins where each element = the number of cm in the [x,y,z] dimension
pbins_3d = 0; % number of bins to pad around data, can also be a 3 element vector i.e. [2 2 4] for anisotropic padding where each element = the number of bins in the [x,y,z] dimension
mdwel_3d = 1; % the minimum dwell time (s) required in each voxel
mdist_3d = 10; % the minimum distince (cm) a voxel has to be from some actual tracking data for it to be calculated
estep_3d = 5; % the step size (pixels) for extending voxels
frcut = 0.2; % percentage cutoff for ratemaps, before detecting place fields
vlcut = 5^3; % the minimum volume (bins) of a place field                    
                
% figure settings
fig_vis = 'off'; % set to 'on' to see figures, 'off' for them to be hidden (this is also faster)
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
[~,~,~] = mkdir([pwd '\klustest3\data']);
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
if ~exist(['klustest3\' cname '_mtint.mat'],'file') || mtint_override
    disp(sprintf('\t...running getDACQDATA'));
    mtint = getDACQDATA(cname,snames);
    disp(sprintf('\t...post-processing mtint'));
    mtint = postprocess_DACQ_data(mtint);
    info = whos('mtint');
    siz = info.bytes / 1000000;
    disp(sprintf('\t...saving %.1fMb mtint',siz))
    save(['klustest3\' cname '_mtint.mat'],'cname','mtint','-v7.3');
    loaded_mtint = 1;
else
    disp(sprintf('\t...loading saved mtint'));
    load(['klustest3\' cname '_mtint.mat'],'-mat');
    loaded_mtint = 1;
end % if ~exist(['klustest3\' cname '_mtint.mat'],'file') || mtint_override
disp(sprintf('\t...done'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get some variables that will be important later
disp(sprintf('Extracting initial data...'))
disp(sprintf('\t...recording sessions: %d',numel(snames)))

%% Get the position data for the sessions
% manually define some settings
samp_rate_hz = 50; % synchronised dacqTrack data should be 50Hz
pos_tb = samp_rate_hz / 1000;    
pixel_ratio = ratio_3d;
leds = 1; % typically we only use 1 LED

pox = []; poy = []; poz = []; % initialise variables
for ss = 1:nsess % for each session
    snow = snames{ss}; % get its filename
    fcheck = [pwd '\triangulation\' snow '_merge_data.mat']; % load its 3d reconstructed path
    if exist(fcheck,'file') % check it exists
        load(fcheck);
    else
        error('Predetermined 3D trajectory not found: %s\n ...exiting',fcheck)
    end % if exist(fcheck,'file')

    posn = datout.weight_mean_merge_smoothed; % extract the smoothed, weighted mean trajectory
    scount_now = mtint.pos.trial_samples(ss); % only take the first n samples, where n = the number of samples in the dacqUSB data (often there is one extra data point in the synchronised data)
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
sess_data = partSESS(sess_data,mtint);
sess_parts = numel(sess_data.names);
disp(sprintf('\t...done'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve and analyse data
disp(sprintf('Analysing tetrode data...'))
%% For every available tetrode
for e = 1:length(tetrodes) 
    tet = tetrodes(e); % tet = the current tetrode
    tet_s = ['t_' num2str(tet)]; % string for use in structure array
    disp(sprintf('\tLoading electrode %d...',tet));
    
    pos_samps = ceil(mtint.tetrode(tet).ts * 50);    
    spik_count = numel(mtint.tetrode(tet).nspike_cut); % retrieve the number of spikes recorded on this tetrode

    disp(sprintf('\t\t...%d spikes detected',spik_count));
    if clusters == 0
        clus = unique(mtint.tetrode(tet).cut); % find the clusters logged for this tetrode
        clus(clus == 0) = []; % find and delete cluster 0 
    else
        clus = clusters;
    end % if clusters == 0
    clus_count = numel(clus);   
    disp(sprintf('\t\t\t...%d data clusters detected',clus_count));
    
%% If there are clusters    
    if clus_count > 0
        disp(sprintf('\t\t\t\t...starting analysis'));
        disp(sprintf('\t\t\t\t...progress: 0%%'));

        % get the spike times for this tetrode
		spiketime = mtint.tetrode(tet).ts;
        
        % calculate the quality of clusters in this session
        clear isod lratio nd cd fdata nfets
        [isods,lratios,~,~,~,~] = clusterQUALITY(cname,tet,clus);        
        
%% For every detected cluster           
        for cc = 1:length(clus)
            clu = clus(cc); % clu = the current cluster
            clu_s = ['c_' num2str(clu)]; % string for use in structure array
            
            isod = isods(cc); % retrieve cluster quality info
            lratio = lratios(cc);
            sdata.(tet_s).(clu_s).cluster_isolation = isod; % add data to structure
            sdata.(tet_s).(clu_s).cluster_lratio = lratio; % add data to structure
           
            % get some vectors that we can use to sort data
			clu_identity = mtint.tetrode(tet).cut; % clu_assign is a vector of numbers, one for each spike, each number corresponds to a cluster
            clu_indx = (clu_identity == clu); % logical vector, one number per spike, 1 if in cluster, 0 if not
			n_spikes = length(find(clu_identity == clu));
			pos_identity = pos_samps; % pos_identity is a vector of numbers, one for each spike, each number corresponds to a position data point
            
            % sort out spike data
            spt = spiketime(clu_indx); % the time point for every spike in this cluster        
            spkindx = knnsearch(pot,spt); % get an index of position data for each spike (i.e. which pox,poy,poz value most closely corresponds to each spike)    
            spx = pox(spkindx);
            spy = poy(spkindx);                       
            spz = poz(spkindx);    
            spx = spx(:);
            spy = spy(:);
            spz = spz(:);
            
            sdata.(tet_s).(clu_s).spx = single(spx); % add data to structure
            sdata.(tet_s).(clu_s).spy = single(spy); % add data to structure
            sdata.(tet_s).(clu_s).spz = single(spz); % add data to structure            
            sdata.(tet_s).(clu_s).spt = single(spt); % add data to structure                      
            hd = double(head_direct(pos_identity(clu_identity == clu))); % get an index of which spikes belong to this cluster, then get an index of whcih positions these spikes correspond to, then use this to get the correct head directions       

            sess_pos = cell(1,sess_parts); % will contain all position data, seperated into goals with an index of trajectories
            sess_spk = cell(1,sess_parts); % will contain all spike data, seperated into goals with an index of trajectories     

%% For each session       
            for ss = 1:sess_parts % for every session named in sess_data
                session_now = sess_data.names{ss}; % the name of the current session
                sessindx = find(sess_data.assignment == ss); % find which time values correspond to this session (maybe non-sequential recordings should be combined)
                session_times = sess_data.times(sessindx,:);
                
                sspk = []; spos = [];            
                session_duration = 0;             
                for rr = 1:length(session_times(:,1)) % for every pair of time values associated with this goal  
                    rec_time = session_times(rr,:); % get the start and end time
                    sindx = find(spt > rec_time(1) & spt < rec_time(2)); % find the spikes falling into this session
                    pindx = find(pot > rec_time(1) & pot < rec_time(2)); % find the position data falling into this session
                                      
                    % concatenate all the data
                    sspk = [sspk; spx(sindx) spy(sindx) spz(sindx) spt(sindx) hd(sindx) ones(size(spx(sindx))).*traj_now];
                    spos = [spos; pox(pindx) poy(pindx) poz(pindx) pot(pindx) ones(size(pox(pindx))).*traj_now];                
                    session_duration = session_duration + (trial_time(2)-trial_time(1));
                end % for tt = 1:length(goal_times(:,1))                                   
                % assign data to arrays
                pox = spos(:,1);
                poy = spos(:,2);
                poz = spos(:,3);                
                pot = spos(:,4);
                spx = sspk(:,1);
                spy = sspk(:,2);
                spz = sspk(:,3);                
                spt = sspk(:,4);    
                sdata.(session_now).pos = spos;

                lname = ['triangulation\' session_now '_lattice.mat'];
                if exist(lname,'file')
                    load(lname);
                else
                    if sess_data.dimensions(ss)
                    [lx,ly,lz] = makeLATTICE_v2(pox,poy,poz,pratio,lsize,lcomp);
                    save(lname,'lx','ly','lz')
                end % if exist('3D matfiles\lattice.mat')
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                if part_dims(ss) == 3
                    if exist('triangulation\lattice.mat','file')
                        load('triangulation\lattice.mat')
                    else
                        [lx,ly,lz] = makeLATTICE(gpox,gpoy,gpoz);
                        save('triangulation\lattice.mat','lx','ly','lz')
                    end % if exist('3D matfiles\lattice.mat')
                    make_lat = 0;
                    
                    % remove extraneous data points
                    offset = 100;
                    nindx = gpox < min(lx)-offset | gpox > max(lx)+offset | gpoy < min(ly)-offset | gpoy > max(ly)+offset | gpoz < min(lz)-offset | gpoz > max(lz)+offset;
                    gpox(nindx) = [];
                    gpoy(nindx) = [];
                    gpoz(nindx) = [];
                    gpot(nindx) = [];
                    gtdx(nindx) = [];
                    nindx = gspx < min(lx)-offset | gspx > max(lx)+offset | gspy < min(ly)-offset | gspy > max(ly)+offset | gspz < min(lz)-offset | gspz > max(lz)+offset;
                    gspx(nindx) = [];
                    gspy(nindx) = [];
                    gspz(nindx) = [];
                    gspt(nindx) = [];
                    ghd(nindx) = [];
                else
                    lx = 0; % dummy data
                    ly = 0;
                    lz = 0; 
                end % if part_dims(gg) == 3
                lx1 = lx;
                ly1 = ly;
                lz1 = lz;
                sdata.lattice_pos = [lx1 ly1 lz1];
                              
                % work out the cluster's goal firing rate
                gfrate = numel(gspx) / session_duration;
                sdata.(tet_s).(clu_s).(session_now).frate = gfrate;
                sdata.(session_now).duration = session_duration;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate 3D ratemap half sessions and correlate them


% keyboard
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                if numel(gspx) > 0 % if there is at least one spike                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate 3D ratemap and dwellmap
                    datin = struct;
                    datin.pox = gpox;
                    datin.poy = gpoy;
                    datin.poz = gpoz;
                    datin.spx = gspx;
                    datin.spy = gspy;
                    datin.spz = gspz;
                    datin.lx = lx;
                    datin.ly = ly;
                    datin.lz = lz;

                    elookup_name = ['triangulation\' cname '_' session_now '_expansion_lookup.mat'];
                    if exist(elookup_name,'file') && elook_override ~= 1
                        load(elookup_name,'expansion_lookup','setin') % load the expansion lookup table, this will speed things up a lot for large ratemaps, also load the settings used for this lookup, to keep things consistent                        
                        elook_override = 0;
                    else
                        setin = struct;
                        setin.pratio = ratio_3d; % pixels per metre, put 1000 for mm accurate 3D reconstruction
                        setin.srate = samp_rate_hz; % sampling rate of data (Hz)
                        setin.bsize = bsize_3d; % starting binsize in cm3, can also be a 3 element vector i.e. [10 10 20] for anisotropic bins where each element = the number of cm in the [x,y,z] dimension
                        setin.pbins = pbins_3d; % number of bins to pad around data, can also be a 3 element vector i.e. [2 2 4] for anisotropic padding where each element = the number of bins in the [x,y,z] dimension
                        setin.mindwell = mdwel_3d; % the minimum dwell time (s) required in each voxel
                        setin.mindist = mdist_3d; % the minimum distince (cm) a voxel has to be from some actual tracking data for it to be calculated
                        setin.exstep = estep_3d; % the step size (pixels) for extending voxels
                        setin.gausmoo = 3; % the number of bins over which to smooth position and spike data
                        expansion_lookup = [];
                        elook_override = 0;
                        disp(sprintf('\b (computing expansion)'));
                    end % if exist(elookup_name,'file') && elook_override ~= 1                    
                    
                    [ratemap,dwellmap,~,~,~,~,expansion_lookup] = generate3Dmap(datin,setin,expansion_lookup);
                    sdata.(tet_s).(clu_s).(session_now).ratemapA = uint16(ratemap);
                    sdata.(session_now).dwellmapA = uint16(dwellmap);                       
                    sdata.(session_now).expansion_map = expansion_lookup;                        
                    save(elookup_name,'expansion_lookup','setin') % save the expansion lookup table, this will speed things up a lot for large ratemaps, also save the settings used for this lookup, to keep things consistent
                    lx2 = expansion_lookup.lx;
                    ly2 = expansion_lookup.ly;
                    lz2 = expansion_lookup.lz;   
                    sdata.lattice_rmap = [lx2 ly2 lz2];
                    thresh = nanmean(nanmean(nanmean(ratemap)));
                    ratemap(ratemap < thresh) = NaN;                         
                    %thresh = 0.2;
                    %ratemap(ratemap < thresh * max2(ratemap)) = NaN;   
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%% Initial ratemap analyses       
                    % spatial information content and sparsity
                    spati = Skaggs3D(ratemap,dwellmap);
                    spars = Spars3D(ratemap,dwellmap);                  
                    sdata.(tet_s).(clu_s).(session_now).sparsity = spars;
                    sdata.(tet_s).(clu_s).(session_now).spatinfo = spati;
                    
                    % find place fields
                    [fields,fdata] = fields3D(ratemap,frcut,vlcut,bsize_3d);                      
                    sdata.(tet_s).(clu_s).(session_now).nfields = fields;
                    sdata.(tet_s).(clu_s).(session_now).field_data = fdata;

                    % skaggs compression matrices
                    [skag1,skag2,skag3,mmap1,mmap2,mmap3] = SkaggsCompress3D(ratemap,dwellmap);
                    sdata.(tet_s).(clu_s).(session_now).spatinfo_dims = [skag1,skag2,skag3];
                    
                    % compression matrix grid scores
                    auto1 = xPearson(mmap1);
                    auto2 = xPearson(mmap2);
                    auto3 = xPearson(mmap3);
                    [gscor1,gspac1,fs1,go1,ge1] = GridAnalysis(auto1,setin.bsize);      
                    [gscor2,gspac2,fs2,go2,ge2] = GridAnalysis(auto2,setin.bsize);      
                    [gscor3,gspac3,fs3,go3,ge3] = GridAnalysis(auto3,setin.bsize);  
                    sdata.(tet_s).(clu_s).(session_now).gscore_comp = [gscor1,gscor2,gscor3]; 
                    sdata.(tet_s).(clu_s).(session_now).gspacing_comp = [gspac1,gspac2,gspac3]; 
                    sdata.(tet_s).(clu_s).(session_now).gfsize_comp = [fs1,fs2,fs3]; 
                    sdata.(tet_s).(clu_s).(session_now).gorientation_comp = [go1,go2,go3]; 
                    sdata.(tet_s).(clu_s).(session_now).gellipticity_comp = [ge1,ge2,ge3]; 
                    
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
                    %[gx,gy,gz,gxm,gym,gzm,gxa,gya,gza] = GridAutoSlice3D(ratemap,setin.bsize,xs,ys,zs);
                    %[sx,sy,sz,sxm,sym,szm,sxd,syd,szd] = SkaggsAutoSlice3D(ratemap,dwellmap,xs,ys,zs);
                    [gx,gy,gz,gxm,gym,gzm,gxa,gya,gza] = GridAutoSlice3D(ratemap,setin.bsize);
                    [sx,sy,sz,sxm,sym,szm,sxd,syd,szd] = SkaggsAutoSlice3D(ratemap,dwellmap);                    
                    sdata.(tet_s).(clu_s).(session_now).gscore_slice = {gx,gy,gz}; 
                    sdata.(tet_s).(clu_s).(session_now).spatinfo_slice = {sx,sy,sz}; 

                    % spike distributions in lattice maze
                    pos = [gpox gpoy gpoz gpot];                    
                    spk = [gspx gspy gspz gspt];                  
                    sdists = spikeDISTRIBUTION(spk,pos,{lx ly lz},0);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1, ratemap, spike plot etc
                        fig_cell = figure('visible',fig_vis,'Position',[100, 100, 1024, 800]);
                        clf % clear it for new data (quicker than opening a new figure each time)
                        set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
                        set(gcf,'color','w'); % makes the background colour white
                        colormap(jet(256)); % to make sure the colormap is not the horrible default one
                        fig_hor = 6; % how many plots wide should it be
                        fig_ver = 6; % how many plots tall should it be
                        fspac = 0.01; % the spacing around the plots, on all sides
                        fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
                        fmarg = 0.03; % the margins around the plots, at the edge of the figure
                        fsiz = 10; % the fontsize for different texts
                        fsiz2 = 10;
                        flw = 2; % the line width for different plots

                        %% add an annotation to the figure with some important info
                        ann_str = sprintf('Tetrode: %d, Cluster: %d, Part: %s, Spikes: %d, Time: %d, Frate: %.2f, Analysed: %s',tet,clu,session_now,numel(gspx),session_duration,gfrate,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
                        annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');                     

% spikes and position plot                        
                        ax1 = subaxis(fig_ver,fig_hor,[1 2 7 8],'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                        plot3(gpox,gpoy,gpoz,'k')
                        hold on
                        plot3(gspx,gspy,gspz,'r.','MarkerSize',20)
                        h = line(lx,ly,lz);
                        set(h,'Color',[0.5 0.5 0.5 0.5],'LineWidth',1,'LineStyle','-')
                        axis off
                        view(3)
                        daspect([1 1 1])
                        xlabel('X');
                        ylabel('Y');
                        zlabel('Z');
                        axis vis3d
                        camproj perspective
                        rotate3d on
                        hold on
                        h = line(lx1,ly1,lz1);
                        set(h,'Color',[0.5 0.5 0.5 0.5],'LineWidth',1,'LineStyle','-');

                        % add annotation with spatial information etc
                        anns = sprintf('Spk:%d\nt:%ds\nFr:%.2f',numel(gspx),round(session_duration),gfrate);                  
                        annotation('textbox',[0.01, 0.9, 1, 0],'string',anns,'FontSize',fsiz,'LineStyle','none','interpreter','none'); 
                        
% ratemap (adaptive binned with cutoff)                     
                        ax2 = subaxis(fig_ver,fig_hor,[3 4 9 10],'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                        vol3d('cdata',ratemap,'texture','3D');
                        alphamap('rampup');
                        %alphamap(0.5 .* alphamap);
                        colormap('jet');
                        caxis([0 nanmax(nanmax(nanmax(ratemap)))]);
                        axis on
                        view(3)
                        daspect([1 1 1])
                        xlabel('X');
                        ylabel('Y');
                        zlabel('Z');
                        axis vis3d
                        camproj perspective
                        rotate3d on
                        title(sprintf('spars=%.1f,skaggs=%.2f',spars,spati),'FontSize',fsiz2)
                        hold on
                        h = line(lx2,ly2,lz2);
                        set(h,'Color',[0.5 0.5 0.5 0.5],'LineWidth',1,'LineStyle','-');

                        % add annotation with spatial information etc
                        anns = sprintf('Mx:%.1f\nMn:%.1f\nSI:%.2f\nSy:%.2f',nanmax2(ratemap),nanmean2(ratemap),spati,spars);                  
                        annotation('textbox',[0.31, 0.9, 1, 0],'string',anns,'FontSize',fsiz,'LineStyle','none','interpreter','none');                          
 
% place field convexhulls                        
                        subaxis(fig_ver,fig_hor,[5 6 11 12],'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                        % plot overall position data outline
                        dwellmap2 = dwellmap;
                        dwellmap2(dwellmap < 1) = NaN; % remove voxels with a coverage less than 1 second (the minimum when generating ratemap)
                        [p2,p1,p3] = ind2sub(size(dwellmap2),find(~isnan(dwellmap2))); % find the subscript indices of these elements
                        K = convhull(p1,p2,p3); % find the convexhull of these points
                        tri = trisurf(K,p1,p2,p3); % plot it
                        set(tri,'FaceColor','k','EdgeAlpha',0.1)
                        alpha(.1)
                        axis on
                        view(3)
                        daspect([1 1 1])
                        xlabel('X');
                        ylabel('Y');
                        zlabel('Z');
                        axis vis3d
                        camproj perspective
                        rotate3d on
                        hold on

                        % plot place field convexhulls
                        annotation('textbox',[0.93, 0.98, 1, 0],'string',sprintf('fields:%d',fields),'FontSize',fsiz,'LineStyle','none','interpreter','none'); 
                        ypos = 0.96;                        
                        
                        for ff = 1:fields % for every detected place field
                            % plot outline of this place field
                            ps = fdata(ff).voxel_index;
                            K = convhull(ps(:,2),ps(:,1),ps(:,3));
                            tri = trisurf(K,ps(:,2),ps(:,1),ps(:,3));
                            set(tri,'FaceColor','r');
                            alpha(.1);

%                             % plot field eigen vectors
%                             eve = fdata(ff).eigen; % eigen vectors
%                             troid = fdata(ff).centroid_weighted; % centroid
%                             axlengths = fdata(ff).axlengths.c; % axis lengths  
%                             plot3([troid(2) troid(2)+eve(2,1)*4],[troid(1) troid(1)+eve(1,1)*4],[troid(3) troid(3)+eve(3,1)*4],'r','LineWidth',3);
%                             plot3([troid(2) troid(2)+eve(2,2)*4],[troid(1) troid(1)+eve(1,2)*4],[troid(3) troid(3)+eve(3,2)*4],'g','LineWidth',3);
%                             plot3([troid(2) troid(2)+eve(2,3)*4],[troid(1) troid(1)+eve(1,3)*4],[troid(3) troid(3)+eve(3,3)*4],'b','LineWidth',3);    
                            
                            % add annotation to the plot showing the skaggs SI content of each plane
                            anns = sprintf('%.2fcm3\necc:%.2f\nelg:%.2f',fdata(ff).volume_cm3,fdata(ff).eccentricity,fdata(ff).elongation1);    
                            annotation('textbox',[0.93, ypos, 1, 0],'string',anns,'FontSize',fsiz,'LineStyle','none','interpreter','none');                              
                            ypos = ypos - 0.07;
                        end % for ff = 1:fields % for every detected place field                               
                        
% skaggs compression matrices
                        subaxis(fig_ver,fig_hor,[13 14 19 20],'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                        % get the sizes of the autocorrelation maps
                        mmap1b = mmap2;
                        mmap2b = fliplr(mmap1);
                        mmap3b = rot90(mmap3,3);
                        [s1a,s2a] = size(mmap1b); % back right wall
                        [s1b,s2b] = size(mmap2b); % back left wall                
                        [s1c,s2c] = size(mmap3b); % floor       
                        
                        % add outlines of a cube
                        xpatch = ([0 1 1 0 0 0; 1 1 0 0 1 1; 1 1 0 0 1 1; 0 1 1 0 0 0]*(s1c-1))+1;
                        ypatch = ([0 0 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 0 1 1 1 1]*(s2c-1))+1;
                        zpatch = ([0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 1 1 0 1; 1 1 1 1 0 1]*(s1a-1))+1;
                        for vv = 1:6
                            h = patch(xpatch(:,vv),ypatch(:,vv),zpatch(:,vv),'w');
                            set(h,'edgecolor','k','FaceColor','none')
                        end % for vv = 1:6     
                        hold on                        
                        
                        % plot the compressed maps onto the sides of the cube                        
                        surf([s1c s1c; s2a s2a],[s2c 1; s2c 1],[s1a s1a; 1 1],'CData',mmap1b,'FaceColor','texturemap'); % back right wall
                        surf([1 s1c; 1 s1c],[s2c s2c; s2c s2c],[s1b s1b; 1 1],'CData',mmap2b,'FaceColor','texturemap'); % back left wall    
                        surf([1 s1c; 1 s1c],[s2c s2c; 1 1],[1 1; 1 1],'CData',mmap3b,'FaceColor','texturemap'); % floor   
                        axis on
                        view(3)
                        daspect([1 1 1])
                        xlabel('X');
                        ylabel('Y');
                        zlabel('Z');
                        axis vis3d
                        camproj perspective
                        rotate3d on
                        
                        % add annotation to the plot showing the skaggs SI content of each plane
                        anns = sprintf('SI\nYZ:%.2f\nXZ:%.2f\nXY:%.2f',skag1,skag2,skag3);                  
                        annotation('textbox',[0.31, 0.6, 1, 0],'string',anns,'FontSize',fsiz,'LineStyle','none','interpreter','none');                     
                        
% compression matrix autocorrelation
                        subaxis(fig_ver,fig_hor,[15 16 21 22],'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                        % get the sizes of the autocorrelation maps
                        auto1b = auto2;
                        auto2b = fliplr(auto1);
                        auto3b = rot90(auto3,3);                        
                        [s1a,s2a] = size(auto1b); % back right wall
                        [s1b,s2b] = size(auto2b); % back left wall                
                        [s1c,s2c] = size(auto3b); % floor    
                        
                        % add outlines of a cube
                        xpatch = ([0 1 1 0 0 0; 1 1 0 0 1 1; 1 1 0 0 1 1; 0 1 1 0 0 0]*(s1c-1))+1;
                        ypatch = ([0 0 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 0 1 1 1 1]*(s2c-1))+1;
                        zpatch = ([0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 1 1 0 1; 1 1 1 1 0 1]*(s1a-1))+1;
                        for vv = 1:6
                            h = patch(xpatch(:,vv),ypatch(:,vv),zpatch(:,vv),'w');
                            set(h,'edgecolor','k','FaceColor','none')
                        end % for vv = 1:6          
                        hold on
                        
                        % plot the autocorrelation maps onto the sides of the cube
                        surf([s1c s1c; s2a s2a],[s2c 1; s2c 1],[s1a s1a; 1 1],'CData',auto1b,'FaceColor','texturemap'); % back right wall
                        surf([1 s1c; 1 s1c],[s2c s2c; s2c s2c],[s1b s1b; 1 1],'CData',auto2b,'FaceColor','texturemap'); % back left wall    
                        surf([1 s1c; 1 s1c],[s2c s2c; 1 1],[1 1; 1 1],'CData',auto3b,'FaceColor','texturemap'); % floor         
                        axis on
                        view(3)
                        daspect([1 1 1])
                        xlabel('X');
                        ylabel('Y');
                        zlabel('Z');
                        axis vis3d
                        camproj perspective
                        rotate3d on

                        % add annotation to the plot showing the grid score for each plane
                        anns = sprintf('Gscore\nYZ:%.2f\nXZ:%.2f\nXY:%.2f',gscor1,gscor2,gscor3);                  
                        annotation('textbox',[0.31, 0.5, 1, 0],'string',anns,'FontSize',fsiz,'LineStyle','none','interpreter','none');                            

% Firing rate histogram along each dimension
                        sax_pos = [25 31; 26 32; 27 33];
                        natin = {'X','Y','Z'};
                        patin = {gpox gpoy gpoz gpot};
                        lx2 = unique(lx);
                        lx2 = lx2(~isnan(lx2));
                        ly2 = unique(ly);
                        ly2 = ly2(~isnan(ly2));                    
                        lz2 = unique(lz);
                        lz2 = lz2(~isnan(lz2));
                        latin = {lx2,ly2,lz2};
                        for sxa = 1:3 % for every histogram plot
                            subaxis(fig_ver,fig_hor,sax_pos(sxa,:),'Spacing',0.05,'Padding',fpadd,'Margin',0.05);
                            %title(sprintf('Smoothed movements in %s',natin{ff}))
                            barh(sdists.(natin{sxa}).spike_bins,sdists.(natin{sxa}).frate_hist,1,'FaceColor','k')
                            lnow = latin{sxa};
                            xlim([0 max(sdists.(natin{sxa}).frate_hist)])
                            ylim([min(patin{sxa})-20 max(patin{sxa})+20]);
                            v2 = axis;
                            hold on
                            for ll = 1:length(lnow)
                                line([v2(1) v2(2)],[lnow(ll) lnow(ll)],'Color','b')
                            end % for ll = 1:length(lnow)
                            axis(v2)
                            set(gca,'LineWidth',flw,'layer','top','FontSize',8) 
                            box off
                            ylabel(sprintf('%s-coordinate (mm)',natin{sxa}),'FontSize',8)
                            xlabel('Fr (Hz)','FontSize',8)                            
                        end % for sxa = 1:3 % for every histogram plot

% Save the figure        
                        id = ['E' num2str(tet) '_C' num2str(clu) '_' session_now];
                        print(fig_cell,'-dpng','-r150',[pwd '\klustest3\figures\' id '.png'])
                        if save_fig
                            set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')');
                            savefig(fig_cell,[pwd '\klustest3\figures\' id '.fig'],'compact');
                        end % if save_fig
                        close(fig_cell);  
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
%% Figure 2 and 3, ratemap and autocorrelation slices
                        data_order = {gxm, gym, gzm; gxa, gya, gza};
                        data_values = {sx, sy, sz; gx, gy, gz};
                        data_names = {'Skaggs','Grid'};
                        for dd = 1:length(data_order(:,1)) % for each type of data we want to plot in slices
                            xdnow = data_order{dd,1};
                            ydnow = data_order{dd,2};
                            zdnow = data_order{dd,3};
                            xvnow = data_values{dd,1};
                            yvnow = data_values{dd,2};
                            zvnow = data_values{dd,3};

                            fig_slice = figure('visible',fig_vis,'Position',[100, 100, 1024, 800]);
                            set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
                            set(gcf,'color','w'); % makes the background colour white
                            colormap(jet(256)); % to make sure the colormap is not the horrible default one
                            fig_hor = 5; % how many plots wide should it be
                            fig_ver = 3; % how many plots tall should it be
                            fspac = 0.01; % the spacing around the plots, on all sides
                            fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
                            fmarg = 0.03; % the margins around the plots, at the edge of the figure
                            fsiz = 10; % the fontsize for different texts
                            fsiz2 = 10;
                            flw = 1; % the line width for different plots

                            %% add an annotation to the figure with some important info
                            ann_str = sprintf('Tetrode: %d, Cluster: %d, Part: %s, Spikes: %d, Time: %d, Frate: %.2f, Analysed: %s',tet,clu,session_now,numel(gspx),session_duration,gfrate,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
                            annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');                     

% plot a cube to demonstrate x axis slicing
                            ax1 = subaxis(fig_ver,fig_hor,1,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);    
                            xpatch = ([0 1 1 0 0 0; 1 1 0 0 1 1; 1 1 0 0 1 1; 0 1 1 0 0 0]);
                            ypatch = ([0 0 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 0 1 1 1 1]);
                            zpatch = ([0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 1 1 0 1; 1 1 1 1 0 1]);                        
                            for vv = 1:6
                                h = patch(xpatch(:,vv),ypatch(:,vv),zpatch(:,vv),'w');
                                set(h,'edgecolor','k','FaceColor','none')
                            end % for vv = 1:6          
                            hold on
                            mnow = xdnow{round(numel(xdnow)/2)};
                            mnow(isnan(mnow)) = 0;
                            surf([0.5 0.5; 0.5 0.5],[0 1; 0 1],[0 0; 1 1],'CData',mnow,'FaceColor','texturemap'); % back right wall
                            axis off
                            view(3)
                            daspect([1 1 1])
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                            axis vis3d
                            camproj perspective
                            rotate3d on
                            set(findall(gca, 'type', 'text'), 'visible', 'on');

% plot a cube to demonstrate y axis slicing
                            ax1 = subaxis(fig_ver,fig_hor,6,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);    
                            xpatch = ([0 1 1 0 0 0; 1 1 0 0 1 1; 1 1 0 0 1 1; 0 1 1 0 0 0]);
                            ypatch = ([0 0 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 0 1 1 1 1]);
                            zpatch = ([0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 1 1 0 1; 1 1 1 1 0 1]);                        
                            for vv = 1:6
                                h = patch(xpatch(:,vv),ypatch(:,vv),zpatch(:,vv),'w');
                                set(h,'edgecolor','k','FaceColor','none')
                            end % for vv = 1:6          
                            hold on
                            mnow = ydnow{round(numel(ydnow)/2)};
                            mnow(isnan(mnow)) = 0;
                            surf([1 0; 1 0],[0.5 0.5; 0.5 0.5],[0 0; 1 1],'CData',mnow,'FaceColor','texturemap'); % back left wall    
                            axis off
                            view(3)
                            daspect([1 1 1])
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                            axis vis3d
                            camproj perspective
                            rotate3d on
                            set(findall(gca, 'type', 'text'), 'visible', 'on');

% plot a cube to demonstrate z axis slicing
                            ax1 = subaxis(fig_ver,fig_hor,11,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);    
                            xpatch = ([0 1 1 0 0 0; 1 1 0 0 1 1; 1 1 0 0 1 1; 0 1 1 0 0 0]);
                            ypatch = ([0 0 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 0 1 1 1 1]);
                            zpatch = ([0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 1 1 0 1; 1 1 1 1 0 1]);                        
                            for vv = 1:6
                                h = patch(xpatch(:,vv),ypatch(:,vv),zpatch(:,vv),'w');
                                set(h,'edgecolor','k','FaceColor','none')
                            end % for vv = 1:6          
                            hold on
                            mnow = zdnow{round(numel(zdnow)/2)};
                            mnow(isnan(mnow)) = 0;
                            surf([1 0; 1 0],[0 0; 1 1],[0.5 0.5; 0.5 0.5],'CData',mnow,'FaceColor','texturemap'); % floor         
                            axis off
                            view(3)
                            daspect([1 1 1])
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                            axis vis3d
                            camproj perspective
                            rotate3d on                        
                            set(findall(gca, 'type', 'text'), 'visible', 'on');

% plot the slices made along the x-axis
                            ax1 = subaxis(fig_ver,fig_hor,[2 3 4 5],'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);                               
                            XM = cellfun(@(x) padmat(x,1,NaN),xdnow,'UniformOutput',false);                   
                            whole_map = rot90(cell2mat(XM),3); 
                            whole_map = cell2mat(XM);                                
                            im = imagesc(whole_map);
                            set(im,'alphadata',~isnan(whole_map))
                            daspect([1 1 1])
                            axis off                                                                                             

% plot the slices made along the y-axis                        
                            ax2 = subaxis(fig_ver,fig_hor,[7 8 9 10],'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);                               
                            YM = cellfun(@(x) padmat(x,1,NaN),ydnow,'UniformOutput',false);                   
                            whole_map = rot90(cell2mat(YM),3);  
                            whole_map = cell2mat(YM);                              
                            im = imagesc(whole_map);
                            set(im,'alphadata',~isnan(whole_map))
                            daspect([1 1 1])     
                            axis off

% plot the slices made along the z-axis                        
                            ax3 = subaxis(fig_ver,fig_hor,[12 13 14 15],'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);                               
                            ZM = cellfun(@(x) padmat(x,1,NaN),zdnow,'UniformOutput',false);                   
                            whole_map = rot90(cell2mat(ZM),3);    
                            im = imagesc(whole_map);
                            set(im,'alphadata',~isnan(whole_map))
                            daspect([1 1 1])
                            axis off

% add annotations and arrows to the plot
                            annotation('textarrow',[0.3 0.9],[0.9 0.9])
                            annotation('textbox',[0.31, 0.93, 1, 0],'string','Slices along x-axis (y-plane)','FontSize',fsiz,'LineStyle','none','interpreter','none');                                                    
                            annotation('textbox',[0.26, 0.93, 1, 0],'string','x = 0','FontSize',fsiz,'LineStyle','none','interpreter','none');                                                    
                            annotation('textarrow',[0.3 0.9],[0.58 0.58])
                            annotation('textbox',[0.31, 0.61, 1, 0],'string','Slices along y-axis (x-plane)','FontSize',fsiz,'LineStyle','none','interpreter','none');     
                            annotation('textbox',[0.26, 0.61, 1, 0],'string','y = 0','FontSize',fsiz,'LineStyle','none','interpreter','none');                                                    
                            annotation('textarrow',[0.3 0.9],[0.26 0.26])
                            annotation('textbox',[0.31, 0.29, 1, 0],'string','Slices down z-axis (horizontal-plane)','FontSize',fsiz,'LineStyle','none','interpreter','none');                          
                            annotation('textbox',[0.88, 0.29, 1, 0],'string','z = 0','FontSize',fsiz,'LineStyle','none','interpreter','none'); 
%                             xcords = [0.26 0.36 0.46 0.57 0.67 0.78 0.88];
%                             for av = 1:length(xvnow)
%                                 annotation('textbox',[xcords(av),0.74,1,0],'string',num2str(xvnow(av)),'FontSize',fsiz,'LineStyle','none','interpreter','none'); 
%                             end % for av = 1:length(xvnow)
%                             for av = 1:length(yvnow)
%                                 annotation('textbox',[xcords(av),0.44,1,0],'string',num2str(yvnow(av)),'FontSize',fsiz,'LineStyle','none','interpreter','none'); 
%                             end % for av = 1:length(xvnow)                            
%                             for av = 1:length(zvnow)
%                                 annotation('textbox',[xcords(av),0.12,1,0],'string',num2str(zvnow(av)),'FontSize',fsiz,'LineStyle','none','interpreter','none'); 
%                             end % for av = 1:length(xvnow)    
                            
                            %% Save the figure        
                            id = ['E' num2str(tet) '_C' num2str(clu) '_' session_now '_' data_names{dd} '_sliced'];
                            print(fig_slice,'-dpng','-r150',[pwd '\klustest3\figures\' id '.png'])
                            if save_fig
                                set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')');
                                savefig(fig_slice,[pwd '\klustest3\figures\' id '.fig'],'compact');
                            end % if save_fig
                            close(fig_slice);                                                         
                        end % for dd = 1:length(data_order(:,1)) % for each type of data we want to plot in slices
                end % if numel(gspx) > 0 % if there is at least one spike                
            end % for gg = 1:length(part_indx) % for each specified goal
            disp(sprintf('\b %.f%%',cc/length(clus)*100));
        end % for c = 1:length(clusters) % for every detected cluster
    end % if clus_count > 0 % if there are clusters
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
































