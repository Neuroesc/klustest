function KlusterAnalysisTINT(tetrodes,clusters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function utilises many minor functions to both analyse and plot cluster data generated 
%   in Tint. If requested, it will output a figure for each cluster, the cluster space of each
%   tetrode, the cross-correlations of ever cluster on a tetrode and a session data structure: sdata.mat.
%   It will also generate an mtint file (mtint.mat) containing all the tetrode and cluster info.
%   KlusterAnalysisTINT(tetrodes,clusters)
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
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manual settings
close all
fclose('all');
mtint_override = 0; % set this to 1 to force the production of a new mtint file (do this if you reclustered in tint for instance)
ignore_pos = 0; % set to 1 if there is no position data

% map settings
map_padd = 20; % (default 20) the number of pixels to pad the maps with
bin_size = 2.5; % (default 4) bin size in cm for calculating the rate map.
map_sigma = 1.8; % (default 0.5) sigma (gaussian standard deviation) to be used for rate and position map smoothing 
min_dwell = 0.1; % (default 0.1) total number of seconds that rat has to be in a bin for it to count
hd_type = 'density'; % (default 'histogram') enter 'density' for a circular kolmogorov smirnov density estimate plot, enter 'histogram' for the traditional histogram polar plot
hd_bins = 64; % (default 64) the number of bins to use when computing HD plot
wave_type = 'bounded'; % (default 'raw') enter 'bounded' for a mean line and shaded standard deviation area, enter 'raw' for 100 random waveforms and lines for mean and standard deviation
time_bins = 2; % (default 2s) time window over which to compute the spike vs time plot
pmap_type = 'histogram'; % (default 'histogram') enter 'tri' for a delauney triangulated surf plot, enter 'histogram' for a binned map, similar to the ratemap
pmap_bs = bin_size*2; % (default bin_size*2) bin size for phase map, should probably be larger than the ratemap etc

% refractory period settings
tau_r = 2; % length of refractory period in ms
tau_c = 0.75; % lockout period of recording system in ms
                
% figure settings
fig_format = 'png'; % output format of saved images
fig_vis = 'off'; % set to 'on' to see figures, 'off' for them to be hidden (this is also faster)
cross_fig = 0; % plot cross-correlograms for all tetrode clusters
clus_fig = 1; % plot cluster space figures for each tetrode
spike_fig = 1; % plot a figure for each cluster showing various measures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% If no tetrodes are specified assume all of them (these will be checked later)
if ~exist('tetrodes','var') || isempty(tetrodes)
    tetrodes = 1:16;
end % if ~exist('tetrodes','var') || ismepty(tetrodes)

%% If no clusters are specified assume all of them
if ~exist('clusters','var') || isempty(clusters)
    clusters = 0;
end % if ~exist('clusters','var') || ismepty(clusters)

%% Start analysis
tic;
disp('----------------------------------------------------------------------------');
disp(sprintf('Running KlusterAnalysisTINT...'))
sdata = struct; % create an empty structure - together with the mtint file this will hold all of the session data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get session data
disp(sprintf('Identifying sessions...'))
[snames,cname] = getSET(pwd);
sdata.session_names = snames; % add data to structure
sdata.combined_name = cname; % add data to structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create directories								                
disp('Preparing directories...')
[~,~,~] = mkdir('figures');
[~,~,~] = mkdir('data');
disp(sprintf('\t...done'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check the tetrodes
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
if ~exist(['data\' cname '_mtint.mat'],'file') || mtint_override
    disp(sprintf('\t...running getDACQDATA'));
    mtint = getDACQDATA(cname,snames,tetrodes);
    if ~ignore_pos
        disp(sprintf('\t...post-processing mtint'));
        mtint = postprocess_DACQ_data(mtint);
    end % if ~ignore_pos
	save(['data\' cname '_mtint.mat'],'mtint','-v7.3');
else
    disp(sprintf('\t...loading saved data'));
	load(['data\' cname '_mtint.mat'],'-mat');
end % ~exist('Data\all_data.mtint','file')
disp(sprintf('\t...done'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get some variables that will be important later
disp(sprintf('Extracting initial data...'))
disp(sprintf('\t...recording sessions: %d',numel(snames)))

% get the total length of the session
duration = key_value('duration',mtint.pos.header,'num');
disp(sprintf('\t...total session time: %ds',duration))
sdata.session_duration = duration; % add data to structure

% get the position data for the whole session
if ignore_pos
    pot = double(mtint.pos.ts); % extract the time stamp of each position value
    pox = ones(size(pot));
    poy = ones(size(pot));
    head_direct = ones(size(pot));
    com_min_x = 0;
    com_min_y = 0;
    leds = 0;
    poss = 0;
    disp(sprintf('\tWARNING: ignoring position data!'));
else
    position = mtint.pos.xy_pixels;	
    posx = position(:,1); % extract just the x coordinates
    posy = position(:,2); % extract just the y coordinates
    pox = double(posx);
    poy = double(posy);
    pot = double(mtint.pos.ts); % extract the time stamp of each position value

    % adjust the data
    com_min_x = min(pox);
    com_min_y = min(-poy);
    pox = pox - com_min_x;
    poy = -poy - com_min_y;
                
    % get the head direction information for the whole session
    head_direct = double(mtint.pos.dir);
    leds = size(mtint.pos.led_pos,2);
    poss = numel(pox);
end % if ignore_pos
sdata.pox = pox; % add data to structure
sdata.poy = poy; % add data to structure
sdata.pot = pot; % add data to structure
sdata.hd = head_direct; % add data to structure
    
% get the pixel ratio (pixels per meter)
pixel_ratio = key_value('pixels_per_metre',mtint.pos.header,'num');

% get the position data sampling rate (should be 50hz) or 0.05s
samp_rate_hz = key_value('timebase',mtint.pos.header,'num');
samp_rate_hz = samp_rate_hz(1,1);
pos_tb = samp_rate_hz / 1000;

disp(sprintf('\t...positions read: %d',poss));
disp(sprintf('\t...tracking LEDs: %d',leds));
disp(sprintf('\t...pixel ratio: %dppm',pixel_ratio));
disp(sprintf('\t...sample rate: %dHz (%.2fs)',samp_rate_hz,pos_tb));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve and analyse data
disp(sprintf('Analysing lfp data...'))
lfp = mtint.lfp(1).lfp(:,1);
Fs = mtint.lfp(1).Fs(1,1);
lfptime = (0:length(lfp)-1)'/Fs; % make a vector for time
if round(max(lfptime)) ~= round(duration)
    disp(sprintf('WARNING: lfp time %.2f does not match session duration %.2f...',max(lfptime),duration));
end % if max(lfptime) ~= duration

% filter lfp to get theta
cutOffFreq = [4 12]; % lower and upper cutoff frequencies in Hz
[b,a] = butter(4,cutOffFreq/(Fs/2)); % Generate 4th order butterworth filter coefficients
lfpfilt = filtfilt(b,a,lfp); % Apply filter to data using zero-phase filtering

mtint.lfp(1).theta = lfpfilt;
mtint.lfp(1).t = lfptime;
disp(sprintf('\t...done'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% partitioning settings
disp(sprintf('Preparing partitions...'))
part_diginput = 0; % set to 1 to use inputs in .inp file
part_sessions = 1; % set to 1 to seperate input sessions into seperate outputs
part_whole = 1;

% sort out array of start and end times
if part_sessions % if we want to divide the output based on input sessions
    
    dvals = mtint.pos.trial_duration;
    dvals = [0; dvals(:)];
    dvals = cumsum(dvals);
    dvals = [dvals circshift(dvals,-1)];
    dvals(end,2) = sum(dvals(:,1));
    dvals = dvals(1:end-1,:);

    part_times = dvals;
    part_indx = num2cell(1:length(dvals(:,1)));
    for ss = 1:numel(snames)
        part_names{ss} = ['s_' snames{ss}];
    end % for ss = 1:numel(snames)
     
elseif part_diginput % if we want to divide the output based on input sessions
    % sort out vector of goal identifiers
    goal_array = mtint.inps.goal; % load goal values
    part_indx{1} = find(goal_array == 1); % find all these values in goal array
    part_indx{2} = find(goal_array == 2); % find all these values in goal array							                                        
    part_indx{3} = find(goal_array == 3); % find all these values in goal array								                                        
    part_indx{4} = find(goal_array == 4); % find all these values in goal array								                                        
    part_indx{5} = cell2mat(part_indx); % find all these values in goal array
    part_names = {'part1' 'part2' 'part3' 'part4' 'Entire'}; 

    part_times = mtint.inps.tstamps; % timestamp values from file
    part_times = part_times(:); % make sure it is a column
    if mod(numel(part_times),2)
        error('ERROR: there is an odd number of digital inputs... exiting')
    end % if mod(x,2)
    part_times = [part_times(1:2:end,:) part_times(2:2:end)]; % place start and end times side by side

else
    part_indx{1} = 1;
    part_times = [0 sum(mtint.pos.trial_duration)];
    part_names{1} = ['s_' cname]; % append session name with an s_ to avoid structure fields beginning with a number (i.e. an error)
    
end % if part_sessions

% if part_whole
%     part_names = {part_names
%     
%     
%     
%     
% end % if part_whole
disp(sprintf('\t...done'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve and analyse data
disp(sprintf('Analysing tetrode data...'))
tet_start = tic;

cell_counts = zeros(max(tetrodes),3);
for e = 1:length(tetrodes) % for every available tetrode
    tet = tetrodes(e); % tet = the current tetrode
    tet_s = ['t_' num2str(tet)]; % string for use in structure array
    disp(sprintf('\tLoading electrode %d...',tet));
    spik_count = numel(mtint.tetrode(tet).pos_sample); % retrieve the number of spikes recorded on this tetrode
    disp(sprintf('\t\t...%d spikes detected',spik_count));
    clus_count = numel(unique(mtint.tetrode(tet).cut));
    disp(sprintf('\t\t...%d data clusters detected',clus_count));
   
	if clus_count > 0 % if there are clusters
        disp(sprintf('\t\t...generating figures'));
        disp(sprintf('\t\t...progress: 0%%'));

		clusters = unique(mtint.tetrode(tet).cut); % find the clusters logged for this tetrode
		%Cl_0 = find(clusters == 0); % find and delete cluster 0
		%clusters(Cl_0) = []; % find and delete cluster 0

%% get the spike times for this tetrode
		spiketime = mtint.tetrode(tet).ts;
        
        c_count = zeros(1,3); % this array will hold a count of the different cell types
		for cc = 1:length(clusters) % for every detected cluster            
            clu = clusters(cc); % clu = the current cluster
            clu_s = ['c_' num2str(clu)]; % string for use in structure array
            
            clear isod lratio nd cd fdata nfets
            [isod,lratio,nd,cd,fdata,nfets] = clusterQUALITY(cname,tet,clu);
            sdata.(tet_s).(clu_s).cluster_isolation = isod; % add data to structure
            sdata.(tet_s).(clu_s).cluster_lratio = lratio; % add data to structure
           
%% get some vectors that we can use to sort data
			clu_identity = mtint.tetrode(tet).cut; % clu_assign is a vector of numbers, one for each spike, each number corresponds to a cluster
            clu_indx = (clu_identity == clu); % logical vector, one number per spike, 1 if in cluster, 0 if not
			pos_identity = mtint.tetrode(tet).pos_sample; % pos_assign is a vector of numbers, one for each spike, each number corresponds to a position data point
			n_spikes = length(find(clu_identity == clu));
            
%% sort out spike data
            spx = pox(pos_identity(clu_indx)); % the x coordinate of every spike in this cluster
            spy = poy(pos_identity(clu_indx)); % the y coordinate of every spike in this cluster
            spt = spiketime(clu_indx); % the time point for every spike in this cluster
            sdata.(tet_s).(clu_s).spx = spx; % add data to structure
            sdata.(tet_s).(clu_s).spy = spy; % add data to structure
            sdata.(tet_s).(clu_s).spt = spt; % add data to structure
            hd = double(head_direct(pos_identity(clu_identity == clu))); % get an index of which spikes belong to this cluster, then get an indext of whcih positions these spikes correspond to, then use this to get the correct head directions       

            goals_pos = cell(1,length(part_indx)); % will contain all position data, seperated into goals with an index of trajectories
            goals_spk = cell(1,length(part_indx)); % will contain all spike data, seperated into goals with an index of trajectories        
            for gg = 1:length(part_indx) % for each specified goal
                goal_now = part_names{gg}; % the name of the current goal
                goal_indx = part_indx{gg}; % an index of time pairs associated with this goal
                goal_times = part_times(goal_indx,:); % the actual time pairs
                
                gspk = [];
                gpos = [];
                sindax = [];
                part_duration = 0;
                for tt = 1:length(goal_times(:,1)) % for every trial associated with this goal
                    traj_now = goal_indx(tt);
                    trial_time = goal_times(tt,:);
                    sindx = find(spt > trial_time(1) & spt < trial_time(2));
                    pindx = find(pot > trial_time(1) & pot < trial_time(2));
                    tspx = spx(sindx);
                    tspy = spy(sindx);
                    tspt = spt(sindx);
                    thd = hd(sindx);
                    tpox = pox(pindx);
                    tpoy = poy(pindx);
                    tpot = pot(pindx);
                    
                    gspk = [gspk; tspx(:) tspy(:) tspt(:) thd(:) ones(size(tspx(:))).*traj_now];
                    sindax = [sindax; sindx(:)];
                    gpos = [gpos; tpox(:) tpoy(:) tpot(:) ones(size(tpox(:))).*traj_now];
                    part_duration = part_duration + (trial_time(2)-trial_time(1));
                end % for tt = 1:length(goal_times(:,1))
                goals_pos{gg} = gpos;
                goals_spk{gg} = gspk;
                
                % assign data to arrays
                gpox = gpos(:,1);
                gpoy = gpos(:,2);
                gpot = gpos(:,3);
                gspx = gspk(:,1);
                gspy = gspk(:,2);
                gspt = gspk(:,3);                
                gtdx = gpos(:,4);
              
                % work out the cluster's goal firing rate
                gfrate = numel(spx) / part_duration;
                sdata.(tet_s).(clu_s).(goal_now).frate = gfrate;
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process the concatenated goal data 
                % create dwell time map
                map_limits = [min(gpox)-map_padd max(gpox)+map_padd min(gpoy)-map_padd max(gpoy)+map_padd];
                [dwellmap,~] = mapDATA(gpox,gpoy,map_limits,bin_size,pixel_ratio);
                dwellmap = dwellmap .* pos_tb;
                dwellmap(dwellmap < min_dwell) = 0;
                dwellmap = imgaussfilt(dwellmap,map_sigma);
                sdata.(tet_s).(clu_s).(goal_now).spatial_dwellmap = dwellmap; % add data to structure

                if numel(gspx) > 0 % if there is at least one spike
%% Waveforms
                    waves{1} = mtint.tetrode(tet).ch1(sindax,:);
                    waves{2} = mtint.tetrode(tet).ch2(sindax,:);
                    waves{3} = mtint.tetrode(tet).ch3(sindax,:);
                    waves{4} = mtint.tetrode(tet).ch4(sindax,:);

                    mean_wav = {}; % will hold mean waveform data
                    std_wav = {}; % will hold standard deviation of waveform data
                    max_wav = repmat(NaN,1,4); % will hold maximum amps of waveforms
                    min_wav = repmat(NaN,1,4); % will hold minimum amps of waveforms
                    width_wav = repmat(NaN,1,4); % will hold waveform widths
                    for w = 1:4 % for every recording channel
                        wav = double(waves{w});
                        ch = squeeze(mean(wav,1));
                        chs = squeeze(std(wav,1));
                        [maxval,maxt] = max(ch);
                        [postminval postmint] = min(ch(maxt:end));
                        postmint = postmint + maxt - 1;
                        width = postmint - maxt;
                        width = width * (1000/50); 
                        width_wav(w) = width;
                        max_wav(w) = maxval;
                        min_wav(w) = postminval;
                        mean_wav{w} = ch;
                        std_wav{w} = chs;
                    end % for w = 1:4 % for every recording channel

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).wave_widths = width_wav; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).wave_maxs = max_wav; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).wave.mins = min_wav; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).wave_stds = std_wav; % add data to structure

%% Firing rate map and analyses
                    % create spike map
                    [spikemap,~] = mapDATA(spx,spy,map_limits,bin_size,pixel_ratio);
                    spikemap = imgaussfilt(spikemap,map_sigma);

                    % create ratemap
                    ratemap = spikemap ./ dwellmap;
                    ratemap(dwellmap == 0) = NaN;

                    % ratemap analysis
                    skaggs = skaggs_info2(ratemap,dwellmap);
                    spars = sparsity(ratemap,dwellmap);
                    cohe = spatial_coherence(ratemap,dwellmap);

                    % place field analysis
                    [fieldd] = getPFIELDS(ratemap);

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).spatial_ratemap = ratemap; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).spatial_information = skaggs; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).spatial_sparsity = spars; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).spatial_coherance = cohe; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).field_count = length(fieldd.fields(:,1)); % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).field_data = fieldd; % add data to structure

%% Grid autocorrelation and analyses
                    % create autocorrelation
                    automap = GridAutoCorr(ratemap);

                    % autocorrelation analysis
                    [grid_score,grid_spacing,field_size,grid_orientation,grid_ellipticity] = GridAnalysis(automap,bin_size);

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).grid_autocorrelation = automap; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).grid_score = grid_score; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).grid_spacing = grid_spacing; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).grid_field_size = field_size; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).grid_orientation = grid_orientation; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).grid_ellipticity = grid_ellipticity; % add data to structure

%% Head direction analysis
                    ai = linspace(0,2*pi,hd_bins)'; % angles for binning
                    if strcmp(hd_type,'density')
                        hd_s = deg2rad(head_direct); % session head direction in radians
                        hd_c = deg2rad(hd); % cell head direction in radians
                        [hd1] = circ_ksdensity(hd_s,ai,[],0.02); % the session head direction       
                        [hd2] = circ_ksdensity(hd_c,ai,[],0.02); % the cell's head direction
                        hd1 = hd1 .* pos_tb; % convert session HD to time
                        hd3 = hd2 ./ hd1; % calculate HD firing rate

                    elseif strcmp(hd_type,'histogram')
                        hd1 = hist(deg2rad(hd),hd_bins); % the session head direction   
                        hd2 = hist(deg2rad(head_direct),hd_bins); % the cell's head direction
                        % nAll(nAll < sample_rate) = nan; % min dwell cutoff
                        hd1 = hd1 .* pos_tb; % convert session HD to time
                        hd3 = hd2 ./ hd1; % calculate HD firing rate
                        fh = fspecial('average',[1 5]);
                        hd3 = imfilter(hd3,fh,'circular','same');

                    end % if strcmp(hd_type,'density')
                    hd1 = hd1 ./ max(hd1); % normalise session hd
                    hd3 = hd3 ./ max(hd3); % normalise cell hd

                    % head direction analyses
                    rayleigh = circ_r(ai,hd3); % rayleigh vector length
                    mx2 = rad2deg(ai(hd3 == max(hd3))); % preferred angle (location of max frate)

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).hd_frate = hd_c; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).hd_rayleigh = rayleigh; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).hd_maximum = mx2; % add data to structure

%% Spikes vs time 
                    bstvals = (min(gpot):time_bins:max(gpot)); % vector of time points at which we should calculate spike probability
                    [bspikes,~] = histc(spt,bstvals);
                    [bsprobs,xibs] = ksdensity(spt,bstvals);

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).spikes_time_histogram = [bstvals(:),bspikes(:)]; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).spikes_time_ksdensity = [xibs(:),bsprobs(:)]; % add data to structure

%% Spike phase analysis
                    % phase preference plot
                    lfpf = mtint.lfp(1).theta;
                    sdata.theta = lfpf; % add data to structure    
                    t = mtint.lfp.t;

                    % estimate theta phase, based on lfp data
                    [phase_out,~,~] = Phase([t lfpf],spt);
                    spp = phase_out(:,2);

                    % bin the theta phase data
                    ai = -pi:0.1:pi;
                    yi = histc(spp,ai);
                    [phprobs,xiph] = ksdensity([spp; (spp+2*pi)],ai);

                    % example cosine wave
                    swav = cos(ai);
                    swav = ((swav ./ max(swav))+abs(min((swav ./ max(swav))))).*max(yi);

                    phmu = circ_mean(spp);
                    phmud = rad2deg(phmu);

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).phase_mean = phmu; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).spike_phase = spp; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).spike_phase_ideal = swav; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).spike_phase_ksdensity = [phprobs(:),xiph(:)]; % add data to structure

%% Phase map and analyses               
                    if strcmp(pmap_type,'tri')
                        warning('off','MATLAB:delaunay:DupPtsDelaunayWarnId')
                        tri = delaunay(spx,spy);
                        phmap = tri;
                    elseif strcmp(pmap_type,'histogram')
                        phasemap = mapDATA3(spx,spy,spp,map_padd,pmap_bs,pixel_ratio);
                        nindx = find(isnan(phasemap));
                        phasemap(nindx) = 0;
                        phasemap = imgaussfilt(phasemap,map_sigma);
                        phasemap(nindx) = NaN;
                        phasemap = phasemap';
                        phmap = phasemap;
                    end % if strcmp(pmap,'tri')

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).phase_map = phmap; % add data to structure

%% Spike autocorrelation - theta analaysis           
                    [corrdata1,tms1] = spk_acorr2(spt,1,500);

                    % fit a decomposing sine wave to the spike autocorrelogram
                    [thetaR,thetaP,thetaIndx,thetaPowr,thetaLin] = getTHETAfit([tms1' corrdata1]);

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).theta_r = thetaR; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).theta_p_value = thetaP; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).theta_index = thetaIndx; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).theta_power = thetaPowr; % add data to structure

%% Spike autocorrelation - refractory period analysis
                    [corrdata2,tms2] = spk_acorr2(spt,1,50);

                    % Calculate refractory contamination
                    % Navratilova and McNaughton (2016) Grids from bands, or bands from grids? An examination of the effects of single unit contamination on grid cell firing fields
                    % Fee, Mitra, Kleinfeld (1996) Automatic sorting of multiple unit neuronal signals in the presence of anisotropic and non-Gaussian variability
                    half_spike = corrdata2(tms2 >= 0); % take only the positive side of spike autocorrelogram
                    half_time = tms2(tms2 >= 0); % take only the positive side of spike autocorrelogram
                    tau_tot = tau_r - tau_c;
                    Ns = n_spikes; % number of cluster spikes
                    lambda = Ns/duration;  % mean firing rate for cluster 
                    RPV = sum((half_spike ~= 0) & (half_time' <= tau_r)); % get the number of refractory period violations

                    % get Poisson confidence interval on number of expected RPVs
                    conf_int = 95; % percent confidence interval
                    [~,interval] = poissfit(RPV,(100-conf_int)/100); 

                    % convert contamination from number of RPVs to a percentage of spikes
                    RPVs = [RPV interval(1) interval(2)];
                    cont_bounds = NaN(1,3);
                    for i = 1:length(RPVs)
                        RPVnow = RPVs(i);
                        RPVT = 2 * tau_tot * Ns; % total amount of time in which an RPV could occur = the usable refractory period (tau_tot) around each spike (*Ns) and on each side of the spike (*2)
                        RPV_lambda = RPV / RPVT; % rate of RPV occurence
                        p =  RPV_lambda / lambda; % estimate of % contamination of cluster

                        % force p to be a real number in [0 1]
                        if isnan(p)
                            p = 0; 
                        end % if isnan(p)
                        if p > 1
                            p = 1; 
                        end   % if p > 1
                        cont_bounds(i) = p;
                    end % for i = 1:length(RPVs)
                    censored_estimate = (tau_c/1000) * spik_count / duration; % estimate the fraction of spikes not detected because of the system lockout

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).refractory_violations = RPV; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).refractory_contamination = cont_bounds(1); % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).refractory_contamination_95_lower_upper_bounds = cont_bounds(2:3); % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).refractory_censoring = cont_bounds(2:3); % add data to structure

%% Mahalanobis distance, cluster quality analyses
                    [isod,lratio,nd,cd,fdata,nfets] = clusterQUALITY(cname,tet,clu);

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).cluster_isolation = isod; % add data to structure
                    sdata.(tet_s).(clu_s).(goal_now).cluster_lratio = lratio; % add data to structure                    

%% Cell type/identity
                    if ignore_pos
                        frate = sdata.(tet_s).(clu_s).(goal_now).frate;
                        max_wav = sdata.(tet_s).(clu_s).(goal_now).wave_maxs;
                        max_wide = sdata.(tet_s).(clu_s).(goal_now).wave_widths;
                        [v,windx] = max(max_wav);
                        width = max_wide(windx);

                        c_type = [];
                        if width > 250 % putative pyramidal cells
                            if frate > 5
                                c_type = 'high pyramidal';
                            elseif frate < 0.1
                                c_type = 'silent cell';
                            else
                                c_type = 'medium pyramidal';
                            end % if frate > 5
                        elseif width < 250 % putative interneurons
                            c_type = 'Interneuron';
                        end % if width > 250 % putative pyramidal cells 

                    else
                        frate = sdata.(tet_s).(clu_s).(goal_now).frate;
                        max_wav = sdata.(tet_s).(clu_s).(goal_now).wave_maxs;
                        max_wide = sdata.(tet_s).(clu_s).(goal_now).wave_widths;
                        skaggs = sdata.(tet_s).(clu_s).(goal_now).spatial_information;
                        gscore = sdata.(tet_s).(clu_s).(goal_now).grid_score;
                        rvect = sdata.(tet_s).(clu_s).(goal_now).hd_rayleigh;

                        [v,windx] = max(max_wav);
                        width = max_wide(windx);

                        c_type = [];
                        if width > 250 % putative pyramidal cells
                            if frate > 5
                                c_type = 'high pyramidal';
                            elseif frate < 0.1
                                c_type = 'silent cell';
                            else
                                if skaggs > 0.5
                                    c_type = 'place cell';
                                    c_count(1) = c_count(1) + 1;
                                else
                                    c_type = 'medium pyramidal';
                                end % if skaggs > 0.5
                            end % if frate > 5
                        elseif width < 250 % putative interneurons
                            c_type = 'Interneuron';
                        end % if width > 250 % putative pyramidal cells

                        if grid_score > 0.2
                            c_type = [c_type ' & grid cell'];
                            c_count(2) = c_count(2) + 1;
                        end % if grid_score > 0.2

                        if rayleigh > 0.2
                            c_type = [c_type ' & hd cell'];
                            c_count(3) = c_count(3) + 1;
                        end % if grid_score > 0.2     
                    end % if ignore_pos

                    % accumulate data
                    sdata.(tet_s).(clu_s).(goal_now).cell_type = c_type; % add data to structure
                end % if numel(gspx) > 0 % if there is at least one spike
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Goal response figure
analysis_fig = 1;
                if analysis_fig
                    fig_cell = figure('visible',fig_vis,'Position',[100, 100, 1024, 800]);
                    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
                    set(gcf,'color','w'); % makes the background colour white
                    colormap(jet(256)); % to make sure the colormap is not the horrible default one
                    fig_hor = 4; % how many plots wide should it be
                    fig_ver = 4; % how many plots tall should it be
                    fspac = 0.01; % the spacing around the plots, on all sides
                    fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
                    fmarg = 0.03; % the margins around the plots, at the edge of the figure
                    fsiz = 5; % the fontsize for different texts
                    flw = 1; % the line width for different plots

                    %% add an annotation to the figure with some important info
                    ann_str = sprintf('Tetrode: %d, Cluster: %d, Part: %s, Spikes: %d, Time: %d, Frate: %.2f, Analysed: %s',tet,clu,goal_now,numel(gspx),part_duration,gfrate,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
                    annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');     

                    %% these vectors control the subplot positions of each plot
                    cumul_pos = 1; % spikes and path
                    dwell_map = 4; % dwell time map
                    frate_map = 2; % firing rate map
                    auto_corr = 3; % grid autocorrelation
                    head_dirn = 5; % head direction
                    field_map = 6; % place field map
                    phas_lock = 7; % phase precession plot
                    phas_mapi = 8; % phase precession map
                    maha_dist = 9; % mahalanobis distance measure
                    clus_spac = 10; % example cluster space                
                    spik_time = [11 12]; % spikes vs time plot
                    wave_form = 13; % primary waveform plot
                    wave_form_multi = 14; % primary waveform plot
                    auto_thet = 15; % spike autocorrelation - theta
                    auto_refr = 16; % spike autocorrelation - refractory period

%% spike plot with black lines for path and red dots for spikes
                    if ~ignore_pos
                        subaxis(fig_ver,fig_hor,cumul_pos,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                        plot(gpox,gpoy,'k')
                        hold on
                        plot(gspx,gspy,'ro','MarkerFaceColor','r','MarkerSize',2)
                        daspect([1 1 1])
                        axis xy off
                        axis([min(gpox)-15 max(gpox)+15 min(gpoy)-15 max(gpoy)+15]);
                        set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
                        title(sprintf('%d spikes (%.2f Hz)',numel(spx),gfrate));

%% dwell time heatmap
                        subaxis(fig_ver,fig_hor,dwell_map,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                        im = imagesc(dwellmap);
                        set(im,'alphadata',logical(dwellmap));
                        daspect([1 1 1])
                        axis xy off
                        set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
                        title(sprintf('%.2fs (%.2f mins)',part_duration,part_duration/60));  

%% firing rate map
                        if numel(gspx) > 0 % if there is at least one spike
                            subaxis(fig_ver,fig_hor,frate_map,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                            im = imagesc(ratemap);
                            set(im,'alphadata',~isnan(ratemap));
                            title('Ratemap')
                            daspect([1 1 1])
                            axis xy off
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);   
                            colorbar
                            title(sprintf('SI %.2f SP %.2f Cohe %.2f',skaggs,(spars*100),cohe));

%% Place field map    
                            subaxis(fig_ver,fig_hor,field_map,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                            im = imagesc(fieldd.binary_ratemap);
                            set(im,'alphadata',~isnan(ratemap));
                            title(sprintf('Fields: %.f',length(fieldd.fields(:,1))));
                            daspect([1 1 1])
                            axis xy off
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  

%% grid autocorrelation
                            subaxis(fig_ver,fig_hor,auto_corr,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                            im = imagesc(automap);
                            set(im,'alphadata',~isnan(automap));
                            daspect([1 1 1])
                            axis xy off
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
                            title(sprintf('G %.2f S %.2f O %.2f E %.2f',grid_score,grid_spacing,grid_orientation,grid_ellipticity));

%% head direction polar plot
                            subaxis(fig_ver,fig_hor,head_dirn,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                            mmp = mmpolar([ai(:); ai(1)],hd1(:),'k:',[ai(:); ai(1)],hd3(:),'b-','FontSize',fsiz,'Grid','on','RGridVisible','off','RTickVisible','off','TTickDelta',20,'RTickLabelVisible','on','TTickLabelVisible','on');
                            set(mmp,'LineWidth',1.5)
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  
                            title(sprintf('r %.2f max %.2f',rayleigh,mx2));
                            
%% spikes vs time plot
                            subaxis(fig_ver,fig_hor,spik_time,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);  
                            bar(bstvals,bspikes,1,'FaceColor','k')
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  
                            ylabel('Frequency') % label y axis

                            yyaxis right
                            plot(xibs,bsprobs,'b','LineWidth',1)
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  
                            v = axis;
                            axis([0 duration v(3) v(4)]);
                            title('Spikes over time')
                            xlabel('Time (uS)') % label x axis
                            ylabel('Probability') % label y axis

%% Phase preference plot           
                            subaxis(fig_ver,fig_hor,phas_lock,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);  
                            yyaxis left
                            ai = -pi:0.1:3*pi;
                            yi = [yi' yi'];
                            ai = ai(:);
                            yi = yi(:);
                            bar(ai,yi);
                            hold on
                            plot(ai,[swav(:); swav(:)],'b:');
                            ylabel('Frequency') % label y axis
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  

                            yyaxis right
                            plot(xiph,phprobs,'k')
                            title('Theta phase')
                            ylabel('Probability') % label y axis
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  

                            xlabel('Theta Phase') % label x axis
                            set(gca,'Xlim',[0,2*pi]);
                            v = axis;
                            set(gca,'Ylim',[0,v(4)]);
                            ax = gca;
                            ax.XTick = linspace(0,2*pi,7); % change Xtick locations to these values
                            ax.XTickLabel = {'0','60','120','180','240','300','360'}; % change Xtick labels to these values
                            title(sprintf('mu: %.2f (%.f)',phmu,phmud));

%% Phase map                
                            subaxis(fig_ver,fig_hor,phas_mapi,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                            if strcmp(pmap_type,'tri')
                                trisurf(phmap,spx,spy,spp);
                                shading('interp');
                                daspect([1 1 1])
                                view(90,90);
                                axis xy off
                                axis([min(gpox)-15 max(gpox)+15 min(gpoy)-15 max(gpoy)+15]);
                                set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);                 
                                title('Phase map')
                                caxis([-pi,pi]);

                            elseif strcmp(pmap_type,'histogram')
                                im = imagesc(phmap);
                                set(im,'alphadata',~isnan(phmap));
                                title('Phasemap')
                                daspect([1 1 1])
                                axis xy off
                                set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);   
                                caxis([-pi,pi]);
                                colorbar

                            end % if strcmp(pmap,'tri')

%% main waveform plot
                            subaxis(fig_ver,fig_hor,wave_form,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');                
                            [~,mindx] = max(max_wav); % find the channel with the highest mean waveform
                            cwaves = double(waves{mindx});

                            rand_index = rand(100,1); % generate a random index that we can use to extract some random waveforms
                            rand_index = ceil(rand_index*length(cwaves(:,1)));

                            wavtime = -200:20:780;
                            if strcmp(wave_type,'raw')
                                plot(wavtime,squeeze(cwaves(rand_index,:)),'k','LineWidth',0.5);
                                hold on;
                                plot(wavtime,squeeze(mean(cwaves,1))+squeeze(std(cwaves)),'r--','LineWidth',1);
                                plot(wavtime,squeeze(mean(cwaves,1))-squeeze(std(cwaves)),'r--','LineWidth',1);
                                plot(wavtime,squeeze(mean(cwaves,1)),'r','LineWidth',1);

                            elseif strcmp(wave_type,'bounded')
                                [hl, hp] = boundedline(wavtime,squeeze(mean(cwaves,1)),squeeze(std(cwaves)),'-k');
                                set(hl,'Color','r') % line color
                                set(hl,'LineStyle','-') % line style
                                set(hl,'LineWidth',1) % line width
                                set(hp,'FaceColor','b') % color of area
                                set(hp,'FaceAlpha',1) % transparency of area

                            end % if strcmp(wave_type,'raw')
                            axis([-200 780 min(min(cwaves)) max(max(cwaves))]);

                            v = axis;
                            text(0,v(4)*0.95,sprintf('max %.1fuV, width %.1fuS',max_wav(mindx),width_wav(mindx)),'FontSize',fsiz)
                            axis xy on square
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
                            xlabel('Time (uS)') % label x axis
                            ylabel('Amplitude (uV)') % label y axis

%% four-waveform plot
                            su1 = subaxis(fig_ver,fig_hor,wave_form_multi,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');    
                            hold on
                            for i = 1:4
                                cwaves = double(waves{i});
                                rand_index = rand(100,1); % generate a random index that we can use to extract some random waveforms
                                rand_index = ceil(rand_index*length(cwaves(:,1)));

                                f2diff = double(max([(max(max(waves{1}))-min(min(waves{1}))) (max(max(waves{2}))-min(min(waves{2})))]));

                                wavtime = -200:20:780;
                                if i == 2
                                    wavtime = wavtime + 1000;
                                elseif i == 3
                                    cwaves = cwaves - f2diff;
                                elseif i == 4
                                    wavtime = wavtime + 1000;
                                    cwaves = cwaves - f2diff-20;
                                end % if i == 2

                                if strcmp(wave_type,'raw')
                                    plot(wavtime,squeeze(cwaves(rand_index,:)),'k','LineWidth',0.5);
                                    hold on;
                                    plot(wavtime,squeeze(mean(cwaves,1))+squeeze(std(cwaves)),'r--','LineWidth',1);
                                    plot(wavtime,squeeze(mean(cwaves,1))-squeeze(std(cwaves)),'r--','LineWidth',1);
                                    plot(wavtime,squeeze(mean(cwaves,1)),'r','LineWidth',1);

                                elseif strcmp(wave_type,'bounded')
                                    [hl, hp] = boundedline(wavtime,squeeze(mean(cwaves,1)),squeeze(std(cwaves)),'-k','alpha');
                                    set(hl,'Color','r') % line color
                                    set(hl,'LineStyle','-') % line style
                                    set(hl,'LineWidth',1) % line width
                                    set(hp,'FaceColor','b') % color of area
                                    if i == mindx
                                        set(hp,'FaceAlpha',1) % transparency of area
                                    else
                                        set(hp,'FaceAlpha',0.5) % transparency of area
                                    end % if i == mindx
                                end % if strcmp(wave_type,'raw')
                            end % for i = 1:4
                            p = get(su1,'pos'); % get position of axes
                            set(su1,'pos',[p(1)-0.05 p(2) p(3) p(4)]) % move the axes slightly

                            title('All waveforms')
                            axis xy on square
                            ax = gca;
                            ax.XTick = []; % change Xtick locations to these values
                            ax.YTick = []; % change Xtick locations to these values
                            ax.XTickLabel = {[]}; % change Xtick labels to these values
                            ax.YTickLabel = {[]}; % change Ytick labels to these values
                            axis([-200, 1780, min([min(min(waves{1})),min(min(waves{2}))])-f2diff-20, max([max(max(waves{1})),max(max(waves{2}))])]);
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
                            v = axis;
                            line([790 790],[v(3) v(4)]) % vertical line
                            line([v(1) v(2)],[min([min(min(waves{1})),min(min(waves{2}))])-10 min([min(min(waves{1})),min(min(waves{2}))])-10]) % horizontal line

%% spike autocorrelation - theta plot           
                            su2 = subaxis(fig_ver,fig_hor,auto_thet,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');    
                            bar(tms1,corrdata1,0.9,'k');
                            v1 = axis;
                            axis([tms1(1) tms1(end) v1(3) v1(4)]);
                            xlabel('Time lag (ms)') % label x axis
                            ylabel('Probability') % label y axis
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
                            set(gcf,'visible','off'); % trying to stop the figure popping up
                            hold on
                            plot(tms1,thetaLin,'r','linewidth',2);
                            title(sprintf('r = %.2f, i = %.2f, p = %.2f',thetaR,thetaIndx,thetaPowr));

%% spike autocorrelation - refractory period plot   
                            su3 = subaxis(fig_ver,fig_hor,auto_refr,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');    
                            bar(tms2,corrdata2,0.9,'k');
                            v1 = axis;
                            axis([tms2(1) tms2(end) v1(3) v1(4)]);
                            xlabel('Time lag (ms)') % label x axis
                            ylabel('Probability') % label y axis
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);                 
                            set(gcf,'visible','off'); % trying to stop the figure popping up
                            title(sprintf('RPV %d, RPVc %.2f, Censored %.2f',RPV,cont_bounds(1),censored_estimate));
                            hold on
                            plot([-tau_r; -tau_r],[0 v(4)],'r','LineWidth',1)
                            plot([tau_r; tau_r],[0 v(4)],'r','LineWidth',1) 

%% mahalanobis distance, cluster quality plot
                            subaxis(fig_ver,fig_hor,maha_dist,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');    
                            if ~isempty(cd) && ~isempty(nd)
                                % get maximum x value
                                max_d = ceil(max([max(nd);max(cd)]));
                                if (isnan(max_d)) % to cope with missing channel data
                                    max_d = 1;
                                end % if (isnan(max_d))

                                xi = linspace(0,max_d,10000); % vector of values where we want to estimate ksdensity
                                [vals1,~,~] = ksdensity(cd,xi,'Support',[0 max_d],'Kernel','epanechnikov','Function','pdf','Bandwidth',0.05);
                                [vals2,~,~] = ksdensity(nd,xi,'Support',[0 max_d],'Kernel','epanechnikov','Function','pdf','Bandwidth',0.05);

                                a1 = area(xi,vals1);
                                set(a1,'FaceColor','b')
                                alpha(.5)
                                hold on
                                a2 = area(xi,vals2);
                                set(a2,'FaceColor','k')
                                alpha(.5)

                                set(gca,'XScale','log');
                                set(gca,'Xlim',[0,max_d]);
                                xlabel('Log(Distance)') % label x axis
                                ylabel('Probability') % label y axis
                                set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
                                v = axis;
                                text(0.2,v(4)*0.95,sprintf('IsoD: %.2f, Lratio: %.2f',isod,lratio),'FontSize',fsiz)
                                set(gcf,'visible','off'); % trying to stop the figure popping up
                            end % if ~isempty(cd) && ~isempty(nd)

%% Cluster space plot  
                            subaxis(fig_ver,fig_hor,clus_spac,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');    
                            clusters = unique(mtint.tetrode(tet).cut);
                            clus_count = numel(clusters);
                            clus_cut = mtint.tetrode(tet).cut;
                            load([cname '\kwiktint\' cname '.kk'],'-mat','kkfet','kkset');
                            fetNames = kkfet.names;
                            fetStr = kkfet.string;
                            fs = find(fetStr == 1);

                            [~,dindx] = sort(max_wav,2,'descend'); % find the order of these values
                            mch1 = dindx(1); % the channel with the maximum amplitude
                            mch2 = dindx(2); % the channel with the second maximum amplitude                    
                            plot_features = 1; % plot the first feature used (should be first principle component)

                            d1 = fdata(:,((mch1-1)*nfets)+plot_features(1)); % get this feature data for this channel
                            d2 = fdata(:,((mch2-1)*nfets)+plot_features(1)); % get this feature data for this channel
                            for c_plot = 1:clus_count % for every cluster on this tetrode
                                cnow = clusters(c_plot);
                                cindx = find(clus_cut == cnow);
                                hold on
                                colplot = [0.5 0.5 0.5 0.5];
                                plot(d2(cindx),d1(cindx),'.','MarkerSize',3,'color',colplot); % plot the clusters
                            end % for cc = 1:clus_count
                            cindx = find(clus_cut == clu);
                            plot(d2(cindx),d1(cindx),'.','MarkerSize',3,'color','r'); % re-plot the current cluster in red to make sure it is on top and stands out

                            axis xy on
                            title([num2str(mch1) ' vs ' num2str(mch2) ' (' num2str(clus_count) ' clusters)'])
                            xlabel(fetNames{fs(plot_features(1))}) % label x axis
                            ylabel(fetNames{fs(plot_features(1))}) % label y axis
                            set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);                              
                            
                        end % if numel(gspx) > 0 % if there is at least one spike
                    end % if ~ignore_pos
                        
%% Save the figure        
                    id = [cname '_E' num2str(tet) '_C' num2str(clu) '_goal_' goal_now];
                    set(gcf,'visible','off')
                    print(fig_cell,'-dpng','-r150',['Figures\' id '.png'])
                    close(fig_cell);                
                
                end % if analysis_fig
            end % for gg = 1:length(part_indx) % for each specified goal
            disp(sprintf('\b %.f%%',cc/length(clusters)*100));
            
        end % for c = 1:length(clusters) % for every detected cluster
        
        disp(sprintf('\t\t...detected %d place cells, %d HD cells, %d grid cells',c_count(1),c_count(3),c_count(2)));
        cell_counts(tet,:) = c_count;
        clear c_count;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster cross-correlation figure
        if cross_fig
            if clus_count > 1
                fig_corr = figure('visible',fig_vis,'Position',[100, 100, 1024, 800]);
                set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
                set(gcf,'color','w'); % makes the background colour white
                colormap(jet(256)); % to make sure the colormap is not the horrible default one
                fig_hor = clus_count+1; % how many plots wide should it be
                fig_ver = clus_count+1; % how many plots tall should it be
                fspac = 0.005; % the spacing around the plots, on all sides
                fpadd = 0.005; % the spacing around the plots, on all sides, this takes more space than fspac though
                fmarg = 0.05; % the margins around the plots, at the edge of the figure
                fsiz = 5; % the fontsize for different texts
                flw = 1; % the line width for different plots

                %% add an annotation to the figure with some important info
                ann_str = sprintf('Session: %s, Tetrode: %d, Spikes: %d, Time: %d, Analysed: %s',cname,tet,clu,spik_count,duration,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
                annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');     

                w_width = 50;
                % plot original clusters
                for cc = 1:clus_count
                    clu = clusters(cc); % clu = the current cluster
                    clu_s = ['c_' num2str(clu)]; % string for use in structure array
                    spt = sdata.(tet_s).(clu_s).spike_spt;

                    subaxis(fig_ver,fig_hor,cc+1,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
                    [tisi,sisi] = spikeINTERVALS(spt,w_width);
                    bar(tisi,sisi,1,'k')
                    xlim([-w_width/2,w_width/2])
                    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
                    title(sprintf('C%d',clu),'FontSize',14)

                    subaxis(fig_ver,fig_hor,cc+(clus_count*cc)+1,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
                    bar(tisi,sisi,1,'k')
                    xlim([-w_width/2,w_width/2])
                    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);
                    ylabel(sprintf('C%d',clu),'FontSize',14) % label y axis
                end % for cc = 1:clus_count

                % plot cross-correlations
                pairs = [nchoosek(clusters,2); fliplr(nchoosek(clusters,2))]; % every possible pair of clusters, including auto-correlations
                for pp = 1:length(pairs(:,1))
                    pnow = pairs(pp,:);
                    spt1 = sdata.(tet_s).(['c_' num2str(pnow(1))]).spike_spt;
                    spt2 = sdata.(tet_s).(['c_' num2str(pnow(2))]).spike_spt;
                    [tisi,sisi] = spikeINTERVALS(spt1,w_width,spt2);
                    sindx = sub2ind([clus_count+1,clus_count+1],pnow(2)+2,pnow(1)+2);

                    subaxis(fig_ver,fig_hor,sindx,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
                    bar(tisi,sisi,1,'b');
                    xlim([-w_width/2,w_width/2]);
                    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
                    ax = gca;
                    ax.XTick = []; % change Xtick locations to these values
                    ax.YTick = []; % change Xtick locations to these values
                    axis off
                end % for pp = 1:length(pairs(:,1))
                %% Save the figure    
                id = [cname '_E' num2str(tet) '_cross-correlograms'];
                print(fig_corr,'-dpng','-r300',['Figures\' id '.png'])
                close(fig_corr);
            end % if clus_count > 1
        end % if cross_fig
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster space figure
        if clus_fig
            fig_clu = figure('visible',fig_vis,'Position',[100, 100, 1024, 800]);

            set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
            set(gcf,'color','w'); % makes the background colour white
            colormap(jet(256)); % to make sure the colormap is not the horrible default one
            fig_hor = 4; % how many plots wide should it be
            fig_ver = 2; % how many plots tall should it be
            fspac = 0.03; % the spacing around the plots, on all sides
            fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
            fmarg = 0.03; % the margins around the plots, at the edge of the figure
            fsiz = 5; % the fontsize for different texts
            flw = 1; % the line width for different plots

            %% add an annotation to the figure with some important info
            ann_str = sprintf('Session: %s, Tetrode: %d, Spikes: %d, Time: %d, Analysed: %s',cname,tet,clu,spik_count,duration,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
            annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');     

            %% cluster space
            [~,~,~,~,fdata,nfets,~] = clusterQUALITY(cname,tet); % get the session feature data for this tetrode
            clusters = unique(mtint.tetrode(tet).cut);
            clus_count = numel(clusters);
            clus_cut = mtint.tetrode(tet).cut;
            load([cname '\kwiktint\' cname '.kk'],'-mat','kkfet','kkset');
            fetNames = kkfet.names;
            fetStr = kkfet.string;
            fs = find(fetStr == 1);

            % I was going to add this as a user argument, but the features are not arranged in a convenient way (i.e. this actually means plot the first and second included - which may be
            % any two features. It doesn't mean plot features 1 and 2, which would be PC1 and PC2. In the .fet file the features are just concatenated with no way to tell which is which
            % I might come back to make this option more flexible - the used features are in fetNames and fetStr, the features themselves are in fdata
            plot_features = [1]; 

            annotation('textbox',[0.55, 1, 1, 0],'string',fetNames{fs(plot_features(1))},'FontSize',15,'LineStyle','none','interpreter','none');     
            epairs = nchoosek(1:4,2); % every possible combination of channel pair  
            for pp = 1:length(epairs)           
                subaxis(fig_ver,fig_hor,pp,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
                c1 = epairs(pp,1); % first channel to plot (x axis)
                c2 = epairs(pp,2); % second channel to plot (y axis)
                d1 = fdata(:,((c1-1)*nfets)+plot_features(1));
                d2 = fdata(:,((c2-1)*nfets)+plot_features(1));

                colmap = jet(clus_count);
                linfo = cell(1,clus_count);
                for cc_plot = 1:clus_count
                    cnow = clusters(cc_plot);
                    cindx = find(clus_cut == cnow);
                    hold on
                    plot(d2(cindx),d1(cindx),'.','MarkerSize',2,'color',colmap(cc_plot,:));
                    axis xy on square
                    title([num2str(c1) ' vs ' num2str(c2)]);
                    xlabel(fetNames{fs(plot_features(1))}); % label x axis
                    ylabel(fetNames{fs(plot_features(1))}); % label y axis
                    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
                    linfo{cc_plot} = ['Cluster ' num2str(cnow)];
                end % for cc = 1:clus_count
            end % for pp = 1:length(epairs)   
            warning('off','MATLAB:legend:IgnoringExtraEntries'); % the legend function will want to warn that there are too many plots (lines)
            [legh,objh,~,~] = legend(linfo,'boxoff');
            M = findobj(objh,'type','Line');
            set(M,'MarkerSize',50);
            set(legh,'Position',[0.5 0.03 0.14 0.44],'FontSize',14);

            %% Save the figure        
            id = [cname '_E' num2str(tet) '_cluster_space'];
            print(fig_clu,'-dpng','-r300',['Figures\' id '.png'])
            close(fig_clu);
        end % if clus_fig
    end % if clus_count > 0 % if there are clusters
end % for e = 1:length(tetrodes) % for every available tetrode         

%% Save the session data structure
save(['data\' cname '_sdata.mat'],'sdata','-v7.3'); % save session data
save(['data\' cname '_mtint.mat'],'mtint','-v7.3'); % save mtint file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finish up
toc1 = toc/60;
disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd '\Figures',' &'');">','Figures folder','</a>'])
disp(['Open ','<a href = "matlab: [s,r] = system(tint_location);">','TINT','</a>'])
disp(sprintf('KlusterAnalysisTINT has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
disp(sprintf('Done.'));
disp('-------------------------------------------------------------------------------------');
































