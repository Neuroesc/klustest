function klustest(tetrodes,clusters,rname,cname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function utilises many minor functions to both analyse and plot cluster data generated 
%   in Tint. If requested, it will output a figure for each cluster, the cluster space of each
%   tetrode, the cross-correlations of every cluster on a tetrode and a session data structure (sdata.mat)
%   It will also generate an mtint file (mtint.mat) containing all the tetrode and cluster info.
%   klustest(tetrodes,clusters,rname,cname)
%
%%%%%%%% Inputs
%   tetrodes = (default = 1:16) the tetrodes to run on in a vector format (i.e. [1 2 3 4 5 6] would run on tetrodes 1 to 6)
%   clusters = (default = 0) the clusters to run on, set this to 0 to run on all clusters
%   rname = the rat name/number as a string - this will be used in the sdata structure so its important to give this
%   cname = (optional) function will automatically look for klustakwiked files named 'kwiktint', but you can change this here (i.e if you chose a custom output name in kwiktint)
%
%   part_config = structure assigned at top of function which specifies the output format of the data, e.g. where NNN = the name you want the data to output as
%         part_config.NNN.method = 2; % method of partition: 1 = combine everything, 2 = take a recording session, 3 = use digital inputs
%         part_config.NNN.intervals = 1; % intervals to include: if method = 1 this does nothing, if method = 2 this should specify which recording sessions to include, if method = 3 this should specify which digital input pairs to include
%         part_config.NNN.inputs = []; % starts empty, this will contain keypress-time pairs, used if method = 3
%         part_config.NNN.times = []; % starts empty, this will specify the start and end time(s) in the combined file for this partition
%         part_config.NNN.dimensions = 2; % 2 if 2D session, 3 if 3D - only required for 3D data
%
%         part_config.FDATA.out_name = 'kwiktint'; % the name to use when finding kwiktint files and saving new files, default is 'kwiktint'
%         part_config.FDATA.rat_name = '852'; % any fields under FDATA are saved along with part_config, but are not used for partitioning, the rat number is a useful thing to include here
%         part_config.FDATA.recording_times = []; % starts empty, will contain the start and end times for every recording session
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
%   22/11/16 adapted from KlusterAnalysisTINT
%   30/03/17 started editing the newer version with more emphasis on sdata structure
%   30/03/17 adapted to new partitioning system using a part_config cell array, function can now split based on any inputs
%   31/03/17 started moving figure work into subfunctions to contain them, each uses the sdata structure
%   01/04/17 finished most of the figure subfunctions
%   02/04/17 added kluspartfig, which generates a summary figure containing all the parts, and correlations between them
%   05/04/17 fixed multiple problems in clusterQUALITY, including wrong lratios and reliance on kwiktint files, now using clusterQUALITY_v2
%   13/04/17 added ability to extract parent filenames from .cut file, meaning the function only needs the kwiktint out name
%   19/04/17 fixed head direction plot and analyses
%   20/04/17 changed from using mapDATA to mapDATA_v2, latter allows for maps with same limits and converted position data necessary for overdispersion
%   20/04/17 added overdispersion calculation (my own)
%   05/05/17 added comments, added easier naming of part_config fields
%   10/05/17 fixed bugs with interval times in part_config, added figKEYS to handle interval times, added exceptions to figCLUST and overdispersion to handle overlapping trials
%   12/05/17 fixed bug in getDACQDATA where keypress times were not accumulated in time, fixed bug in partSESS where incorrect interval vector was being loaded
%   12/05/17 added ability to send emails when completed
%   17/05/17 added possibility to calculate shuffled spatial measures
%   19/05/17 replaced read_key with saveKEY, this automatically saves a text file version of the .inp files as it loads them
%   19/05/17 replaced getPARTconfig with getPARTconfig_v2, this saves a much nicer text file which is easier to edit
%   07/06/17 fixed a problem in getDACQDATA where it would try to run on all tetrodes even if a subset is specified
%   22/06/17 fixed error in time allocation
%   10/07/17 changed part_config specification so that it is done in a loop, added information required by contest and klustest3
%   25/07/17 added variable pixel ratio capability (big job!) getDACQDATA now uses getDACQDATAHEADERS instead of the horible key_value functions
%   01/11/17 fixed bug where LFP and session duration did not match
%   01/11/17 replaced celltype with getCELLTYPE, replaced mapDATA_v3 with mapDATA4, replaced GridAnalaysis with gridSCORE2, replaced AutoCorr with ndautoCORR
%   07/11/17
%   v4.1.0
%
% to do: change spike phase to spike times, rather than closest pot time
%
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manual settings
% cname is the name to use when seeking cluster cut data (i.e. the name of your cut file minus the extension).
% specified in the next if section is the value of cname to use
% I have set up kwiktint so it outputs stuff named 'kwiktint' by default in most cases, so this should be the default here. 
% If you changed this in kwiktint you should change it here though.
% If you want cname to have a different value for a one off, pass it in as an input argument or initialise cname here instead:
% cname = 'b1'; % comment/edit/uncomment me when needed
% either way these next lines will be skipped and your input value will be used
if ~exist('cname','var') || isempty(cname) || all(isnan(cname))
    cname = 'kwiktint'; % (default = 'kwiktint') please don't edit me!
end

% rname is the rat name to use when saving data inside the sdata structure
% specified in the next section is the default value of rname to use, just as a filler
% You should pass it in as an input argument or initialise it here with the current rat's name 
% (it's a pain I know, but it's the only way to get the data in):
% rname = '775'; % comment/edit/uncomment me when needed
pname = pwd; sindx=strfind(pname,'\'); rname=pname(sindx(end-1)+1:sindx(end-1)+3); % I can extract the rat name from the file path because I know my folder structure
% either way these next lines will be skipped and your input value will be used
if ~exist('rname','var') || isempty(rname) || all(isnan(rname))
    rname = '000'; % (default = '000') please don't edit me!
end
% rname will automatically have an 'r' appended before it because structures cannot have field names that start with a number

close all
fclose('all');

%% ###################################### %%
% Specify partition configuration settings
part_names = {'square1','lattice','square2'}; % the names you want the outputs to be saved as, these cannot start with a number: i.e. 'session1'
% method of partition, corresponding to each of the names above: 1 = combine everything, 2 = take a recording session, 3 = use digital inputs
part_methods = [2 2 2]; % i.e. if part_methods=[2 2 2 1], we will have 4 parts, the first 3 correspond to some combination of recording sessions (there can be multiple ones) and the last one will include all data
% cell array of vectors indicating which intervals to include in each partition: if method = 1 this does nothing, if method = 2 this should specify which recording sessions to include, if method = 3 this should specify which digital input pairs to include (inf = use all)
part_intervals = {1 2 3}; % i.e. if part_methods=[2 2 2], then if part_intervals={1 2 3}, rec 1 will go in part 1, rec 2 in part 2 and rec 3 in part 3    OR     if part_methods=[2 2 2 2 2], then if part_intervals={1 [2 4 5] 3}, rec 1 will go in part 1, rec 2,4 and 5 in part 2 and rec 3 in part 3
dimensions = [2 3 2]; % for klustest3
% part_etypes = [1 3 2 2 2 2 2 2 2 3 1]; % for contest: environment order: 0 = ignore, 1 = circle, 2 = contextbox, 3 = empty contextbox
% part_rtypes = [0 0 0 90 0 270 0 270 180 0 0]; % for contest: clockwise rotation angle
for pp = 1:length(part_names)
    part_config.(part_names{pp}).method = part_methods(pp);
    part_config.(part_names{pp}).intervals = cell2mat(part_intervals(pp));
    part_config.(part_names{pp}).inputs = []; % starts empty, this will contain keypress-time pairs, used if method = 3
    part_config.(part_names{pp}).times = []; % starts empty, this will specify the start and end time(s) in the combined file for this partition
    part_config.(part_names{pp}).dimensions = dimensions(pp); % for klustest3
%     part_config.(part_names{pp}).etype = part_etypes(pp); % for contest
%     part_config.(part_names{pp}).rtype = part_rtypes(pp); % for contest
end
% any fields under FDATA are saved along with part_config, but are not used for partitioning
part_config.FDATA.out_name = cname; % the name to use when finding kwiktint files and saving new files, default is kwiktint
part_config.FDATA.rat_name = rname; % the rat number is a useful thing to include here as it cannot be retrieved anywhere else
part_config.FDATA.recording_times = []; % starts empty, will contain the start and end times for every recording session
% part_config.FDATA.enames = {'circ','cbox','ecbox'}; % for contest: these environment names correspond to etype
% part_config.FDATA.pratio = [395 446 446 446 446 446 446 446 446 446 395]; % recording session pixel ratios, if this is not empty these will be used instead of the dacq data file values. This may not correspond to your parts unless you are partitioning by method 2
%% ###################################### %%

% overrides
pconfig_override = 0; % set to 1 if you want to ignore and overwrite an existing part_config, set to 2 to run with current part_config settings without overwriting anything
maintain_mtint = 0; % set to 1 to save/load mtint in the base workspace, this saves time when running the function mutliple times (for instance in debugging) but should otherwise be set to 0
mtint_override = 0; % set this to 1 to force the production of a new mtint file (do this if you reclustered in tint for instance)

% map settings
config.rmethod = 'gaussian'; % (default 'nearest') the mapping approach to use, either 'nearest','gaussian','adaptive','KDE'
config.map_padd = 2; % (default 2) the number of bins to pad spatial maps with
config.bin_size = 2; % (default 2) bin size in cm for calculating the rate map (make sure data is in cm or pm/sm values are given)
config.map_sigma = 1.5; % (default 1.5) used by nearest and KDE method, sigma of gaussian to use when smoothing traditional dwell and spike maps, or used as bandwidth of KDE
config.min_dwell = 0.1; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
config.g_sigma = 5; % (default 10) only used by gaussian method - sigma of gaussian used to weight position point and spike point distances, bigger means smoother maps
config.min_dist = 2; % (default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
config.srate = 50; % (default 50) sampling rate of data in Hz, used to calculate time

config.time_bins = 2; % (default 2s) time window over which to compute the spike vs time plot
config.frcut = 0.2; % relative minimum firing rate (% of ratemap max) to be considered a field
config.arcut = 9; % minimum number of pixels to be considered a field, typically people use 9 contiguous pixels
config.minfr = 1; % (Hz) absolute minimum cutoff firing rate to be considered a field
config.over_smooth = 13; % number of position data points over which to smooth instantaneous firing rate when calculating overdispersion
            
% HD settings
config.hd_type = 'density'; % (default 'histogram') enter 'density' for a circular kolmogorov smirnov density estimate plot, enter 'histogram' for the traditional histogram polar plot
config.hd_bins = 64; % (default 64) the number of bins to use when computing HD plot
config.hd_sigma = 0.04; % (default 2) the standard deviation of the gaussian kernel used in HD circular density estimate
config.hd_boxcar = 3; % (defualt 3) the number of bins over which to compute the HD histogram boxcar average

% refractory period settings
tau_r = 2; % length of refractory period in ms
tau_c = 0.75; % lockout period of recording system in ms

% figure settings
run_figCLUST = 1; % set to 1 for individual cluster summary figures
run_figPARTS = 0; % set to 1 for cluster summary figure with ratemaps etc for every part
run_figCROSS = 0; % set to 1 for cross correlogram firing figures
run_figCSPACE = 0; % set to 1 to generate a cluster space plot for each cell
fig_vis = 'off';
save_fig = 0; % set to 1 to save .fig files for 3D viewing

% alert settings
alert_me = 0; % set to 1 for Matlab to email you when the function has completed
alert_email = 'example@gmail.com'; % this is the email address matlab will send the alert to

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start analysis
tic;
disp('----------------------------------------------------------------------------');
disp(sprintf('Running klustest...'))
sdata = struct; % create an empty structure - together with the mtint file this will hold all of the session data
d_start = datestr(now);

% deal with part_config files
[~,~,~] = mkdir([pwd '\klustest\' cname]);
pconfig_name = ['klustest\' cname '\' cname '_part_config.txt'];
if ~exist(pconfig_name,'file') || pconfig_override == 1
    disp(sprintf('\t...writing new part_config'));
    part_config = getPARTconfig_v2(pconfig_name,part_config); % save the current part_config cell array
elseif pconfig_override == 2
    disp(sprintf('\t...writing temp part_config'));
    pconfig_name = ['klustest\' cname '\temp.txt'];
    part_config = getPARTconfig_v2(pconfig_name,part_config); % save the current part_config cell array
else
    disp(sprintf('\t...loading part_config'));
    part_config = getPARTconfig_v2(pconfig_name); % load an existing part_config file
end

% If no output name is specified, assume kwiktint output files with name 'kwiktint'
if isfield(part_config.FDATA,'out_name')
    cname = part_config.FDATA.out_name;
end
sdata.settings.combined_name = cname; % add data to structure

% if no rat name is given use a default number
if isfield(part_config.FDATA,'rat_name')
    rname = num2str(part_config.FDATA.rat_name);
    disp(sprintf('\t...rat name: %s',rname))    
end
if ~exist('rname','var') || isempty(rname) || all(isnan(rname))
    rname = '111';
    disp(sprintf('\t...WARNING: no rat name given, using default: %s',rname))
end
sdata.settings.rat_num = rname;

% If no tetrodes are specified assume all of them
if ~exist('tetrodes','var') || isempty(tetrodes)
    tetrodes = 1:16;
end

% If no clusters are specified assume all of them
if ~exist('clusters','var') || isempty(clusters)
    clusters = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get session data
% Create directories								                
disp('Preparing directories...')
[~,~,~] = mkdir([pwd '\klustest\' sdata.settings.combined_name]);
disp(sprintf('\t...done'))

%% Get session names from cut file
disp(sprintf('Identifying sessions...'))
for tt = 1:length(tetrodes) % for every hypothetical tetrode
    cutname1 = [cname '_' num2str(tetrodes(tt)) '.cut']; % the cut file name
    if any(exist(cutname1,'file')) % see if it exists
        break % if it exists use this filename
    end
end
disp(sprintf('\t...reading %s',cutname1));
[~,etext] = getcut(cutname1); % get data from cutfile, specifically the line which identifies the parent sessions
idx = strfind(etext,': '); % should contain two values
flist = textscan(etext(idx(1)+2:idx(2)-8),'%s','delimiter',','); % get the file name parts of the string
snames = flist{1,1}; % extract the relevant cell
nsess = numel(snames); % session number
disp(sprintf('\t...working on sessions: %s',strjoin(snames,', ')));

% Check tetrodes
disp(sprintf('Assessing data...'))
[tetrodes,mvalue] = getTRODES(snames,tetrodes);
if isempty(mvalue)
    disp(sprintf('\t...tetrodes: %s accounted for',mat2str(tetrodes)))
else
    disp(sprintf('\t...WARNING: tetrodes: %s are incomplete and will be skipped',mat2str(mvalue)))
end
disp(sprintf('\t...done'))

% collect data
sdata.settings.session_names = snames; % add data to structure
sdata.settings.num_sessions = nsess;
sdata.settings.config = config;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read all dacq data or load it
disp(sprintf('Fetching DACQ data...'));
if ~maintain_mtint || mtint_override
    clear mtint
end
if any(strcmp(evalin('base','who'),'mtint')) && ~mtint_override && maintain_mtint
    disp(sprintf('\t...using mtint held in memory'));
    mtint = evalin('base','mtint');
elseif ~exist(['klustest\' cname '\' cname '_mtint.mat'],'file') || mtint_override
    mtint = getDACQDATA(cname,snames,tetrodes); % my replacement for readalldacqdata
    disp(sprintf('\t...post-processing mtint'));
    
    % deal with user defined pixel ratio before we calculate speed etc
    if isfield(part_config.FDATA,'pratio') && ~isempty(part_config.FDATA.pratio)
        hdata = mtint.pos.header;
        if length(part_config.FDATA.pratio) ~= length(hdata)
            disp(sprintf('\tERROR: number of user specified pixel ratios (%d) does not equal number of recordings (%d)... exiting',length(part_config.FDATA.pratio),length(hdata)));
            return
        else
            disp(sprintf('\tWARNING: manual override of pixel ratio values: %s',mat2str(part_config.FDATA.pratio)));            
        end
        
        for pp = 1:length(hdata)
            if isfield(hdata(pp),'num_pos_samples')
                pnow = ones(hdata(pp).num_pos_samples,1) .* part_config.FDATA.pratio(pp);
                hdata(pp).pixels_per_metre_vec = uint16(pnow);
                hdata(pp).pixels_per_metre = part_config.FDATA.pratio(pp);
            end
        end   
        mtint.pos.header = hdata;
    end    

    % post process position data etc
    mtint = postprocessDACQDATA(mtint); % my replacement for postprocess_DACQ_data
    
    % save mtint
    info = whos('mtint');
    siz = info.bytes / 1000000;
    disp(sprintf('\t...saving %.1fMb mtint',siz))
    save(['klustest\' cname '\' cname '_mtint.mat'],'cname','mtint','-v7.3');
else
    disp(sprintf('\t...loading saved mtint'));
    load(['klustest\' cname '\' cname '_mtint.mat'],'-mat');
end
disp(sprintf('\t...done'));
sdata.elecdata.elecdata = mtint.tetrode;
sdata.settings.date = mtint.header_date;
part_config.FDATA.recording_times = mtint.pos.trial_duration;

if maintain_mtint
    assignin('base','mtint',mtint); % leave mtint in base workspace
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% partitioning settings
disp(sprintf('Preparing partitions...'))
% determine partition time values
disp(sprintf('\t...finding part_times'));
part_config = partSESS(part_config,mtint); % process the part_config to get the start and end time of each partition
getPARTconfig_v2(pconfig_name,part_config); % save the current part_config cell array
nparts = length(fieldnames(part_config))-1;
part_names = fieldnames(part_config);
disp(sprintf('\t...done'));
sdata.settings.part_config = part_config;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get some variables that will be important later
disp(sprintf('Extracting initial data...'))
disp(sprintf('\t...recording sessions: %d',numel(snames)))

% get the total length of the session
duration = mtint.pos.total_duration;
disp(sprintf('\t...total session time: %ds',duration))
sdata.settings.session_duration = duration; % add data to structure

% get the position data (dacqUSB data) for the whole session
position = mtint.pos.xy_pixels;	
posx = position(:,1); % extract just the x coordinates
posy = position(:,2); % extract just the y coordinates
pox = double(posx);
poy = double(posy);
pot = double(mtint.pos.ts); % extract the time stamp of each position value
pov = double(mtint.pos.speed(:)); % the running speed throughout the session
poh = double(mtint.pos.dir(:,1)); % HD data
ppm = {mtint.pos.header.pixels_per_metre_vec}.'; ppm = double(vertcat(ppm{:,:})); % pixels per metre

% adjust the data
com_min_x = min(pox);
com_min_y = min(-poy);
pox = pox - com_min_x;
poy = -poy - com_min_y;

% accumulate data
sdata.position.pox = pox;
sdata.position.poy = poy;
sdata.position.pot = pot; 
sdata.position.pov = pov;
sdata.position.poh = poh;
sdata.position.ppm = ppm;
poss1 = numel(pox);
leds = size(mtint.pos.led_pos,2);

% get the position data sampling rate (should be 50hz) or 0.05s
samp_rate_hz = mtint.pos.header(1).sample_rate_num;
samp_rate_hz = samp_rate_hz(1,1);
pos_tb = 1 / samp_rate_hz;
config.srate = samp_rate_hz;

% display results
disp(sprintf('\t...positions read: %d',poss1));
disp(sprintf('\t...tracking LEDs: %d',leds));
disp(sprintf('\t...median pixel ratio: %dppm',nanmedian(ppm)));
disp(sprintf('\t...sample rate: %dHz (%.2fs)',samp_rate_hz,pos_tb));

% accumulate data
sdata.position.position_timebase = pos_tb;
sdata.position.position_srate_hz = samp_rate_hz;
sdata.position.leds = leds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve and analyse data
disp(sprintf('Analysing lfp data...'))
lfp = double(mtint.lfp(1).lfp(:,1));
Fs = mtint.lfp(1).Fs(1,1);
lfptime = (0:length(lfp)-1)'/Fs; % make a vector for time
disp(sprintf('\t...samples read: %d',numel(lfp)));
disp(sprintf('\t...sample rate: %dHz',Fs));
if round(max(lfptime)) ~= round(duration)
    disp(sprintf('\tWARNING: lfp time %.2f does not match session duration %.2f...',max(lfptime),duration));
end

% filter lfp to get theta
disp(sprintf('\t...extracting theta'));
cutOffFreq = [4 12]; % lower and upper cutoff frequencies in Hz
[b,a] = butter(4,cutOffFreq/(Fs/2)); % Generate 4th order butterworth filter coefficients
lfpfilt = filtfilt(b,a,lfp); % Apply filter to data using zero-phase filtering
mtint.lfp(1).theta = lfpfilt;
mtint.lfp(1).t = lfptime;

% get theta phase corresponding to every position sample
lfpf = mtint.lfp(1).theta;
sdata.lfpdatas.theta = int16(lfpf); % add data to structure    
t = mtint.lfp.t;
[phase_out,~,~] = Phase([t lfpf],double(pot));
pop = phase_out(:,2);    
sdata.position.pop = pop;
disp(sprintf('\t...done'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve and analyse data
disp(sprintf('Analysing tetrode data...'))
%% For every available tetrode
for e = 1:length(tetrodes) 
    tet = tetrodes(e); % tet = the current tetrode
    tetstr = ['t' num2str(tet)];
    
    % load this tetrodes data
    disp(sprintf('\tLoading tetrode %d...',tet));
    spik_count = mtint.tetrode(tet).nspike_cut; % retrieve the number of spikes recorded on this tetrode
    spiketime = mtint.tetrode(tet).ts;
    
    % get positions corresponding to spikes
    spikepos = knnsearch(single(pot),single(spiketime)); % the position point for every spike on this tetrode
    sdata.elecdata.elecdata(tet).posindx = spikepos; % index same length as spiketime, numbers correspond to values in pox,poy,pot, can be used to get x,y,t of all spikes
       
    % get a vector of clusters we want to analyse
    if clusters == 0
        clus = unique(mtint.tetrode(tet).cut); % find the clusters logged for this tetrode
    else
        clus = clusters;
    end
    
    % check to see if there are any clusters
    clus_count = numel(clus);   
    disp(sprintf('\t\t...%d spikes detected',spik_count));
    disp(sprintf('\t\t\t...%d data clusters detected',sum(clus~=0))); 
    if ~clus_count || ~any(clus) % if there are no clusters, or if there is only a noise cluster
        continue % skip analysis
    end
    
    % calculate the quality of clusters in this session
    disp(sprintf('\t\t\t\t...starting analysis'));
    disp(sprintf('\t\t\t\t...getting cluster quality'));    
    fetname = [sdata.settings.combined_name '.fet.' num2str(tet)];
    cutname = [sdata.settings.combined_name '_'  num2str(tet) '.cut'];
    fetdata = clusterQUALITY_v2(fetname,cutname);        
    sdata.elecdata.(tetstr).fetdata = fetdata;
    sdata.elecdata.(tetstr).clusters = clus;
    clear fetdata
    
    % get the channel waveforms for this tetrode
    disp(sprintf('\t\t\t\t...getting waveforms'));        
    waves = cell(1,4);
    for ggnow = 1:length(sdata.settings.session_names)
        fnamen = sdata.settings.session_names{ggnow};
        [~,c1,c2,c3,c4] = getspikes([fnamen,'.',num2str(tet)]);   
        waves{1} = [waves{1}; c1];
        waves{2} = [waves{2}; c2];
        waves{3} = [waves{3}; c3];
        waves{4} = [waves{4}; c4];  
    end
         
%% For every detected cluster   
    disp(sprintf('\t\t\t\t...progress: 0%%'));
    for cc = 1:length(clus) % for every cluster
        clu = clus(cc); % clu = the current cluster
        
        % skip noise cluster
        if ~clu % if this is the noise cluster don't continue any further
            disp(sprintf('\b %.f%%',cc/length(clus)*100));
            continue
        end
        
        % prepare unique cell identifier
        uci = ['r' rname '_' sdata.settings.date '_t' num2str(tet) '_c' num2str(clu)]; % unique cell identifier - a string with [rat name, session date, electrode, cluster];
                      
        % get some vectors that we can use to sort data
        clu_identity = mtint.tetrode(tet).cut; % clu_assign is a vector of numbers, one for each spike, each number corresponds to a cluster
        clu_indx = (clu_identity == clu); % logical vector, one number per spike, 1 if in cluster, 0 if not
        cspt = spiketime(clu_indx); % the time point for every spike in this cluster 
        sdata.celldata.(uci).spike_index = logical(clu_indx); % a logical index into all spike time vector, 1 means in this cluster
        
%% For each part       
        for pp = 1:nparts % for every partition
            part_now = part_names{pp}; % the name of the current part
            part_times = part_config.(part_now).times; % the time pairs (intervals) corresponding to this part
                    
            % find what data falls into the intervals associated with this part
            sindax = spiketime > part_times(:,1) & spiketime < part_times(:,2); % logical array, same length as spiketime, 1 means spike is in this part
            sindax = clu_indx & sindax; % logical array, same length as spiketime, 1 when spike is in this cluster and is in this part
            pindax = pot > part_times(:,1) & pot < part_times(:,2); % logical array, same length as pox, 1 means position sample falls into this part
            part_duration = sum(pindax(:))*pos_tb;

            % assign data to sdata
            sdata.position.(part_now).position_index = pindax; % this index tells us which position samples belong in this part
            sdata.position.(part_now).duration = part_duration; % total length of part
            sdata.celldata.(uci).(part_now).part_spike_index = sindax; % this index tells us which spikes belong in this part, this can also be used to find the positions

            % extract position data now - I'm doing this using sdata to demonstrate how it can be done
            ppox = sdata.position.pox ( sdata.position.(part_now).position_index ); % position x for this part
            ppoy = sdata.position.poy ( sdata.position.(part_now).position_index ); % position y for this part
            ppot = sdata.position.pot ( sdata.position.(part_now).position_index ); % position time for this part
            ppov = sdata.position.pov ( sdata.position.(part_now).position_index ); % position velocity for this part
            ppoh = sdata.position.poh ( sdata.position.(part_now).position_index ); % position hd for this part
            pppm = sdata.position.ppm ( sdata.position.(part_now).position_index ); % position pixel ratio for this part
            ppop = sdata.position.pop ( sdata.position.(part_now).position_index ); % position theta phase for this part

            % extract spike data now - I'm doing this using sdata to demonstrate how it can be done
            pspx = sdata.position.pox ( sdata.elecdata.elecdata(tet).posindx ( sdata.celldata.(uci).(part_now).part_spike_index )); % spike x for this part
            pspy = sdata.position.poy ( sdata.elecdata.elecdata(tet).posindx ( sdata.celldata.(uci).(part_now).part_spike_index )); % spike y for this part
            pspt = sdata.position.pot ( sdata.elecdata.elecdata(tet).posindx ( sdata.celldata.(uci).(part_now).part_spike_index )); % spike t for this part
            pspt = spiketime(sindax);
            pspv = sdata.position.pov ( sdata.elecdata.elecdata(tet).posindx ( sdata.celldata.(uci).(part_now).part_spike_index )); % spike velocity for this part            
            psph = sdata.position.poh ( sdata.elecdata.elecdata(tet).posindx ( sdata.celldata.(uci).(part_now).part_spike_index )); % spike hd for this part
            pspm = sdata.position.ppm ( sdata.elecdata.elecdata(tet).posindx ( sdata.celldata.(uci).(part_now).part_spike_index )); % spike pixel ratio for this part
            pspp = sdata.position.pop ( sdata.elecdata.elecdata(tet).posindx ( sdata.celldata.(uci).(part_now).part_spike_index )); % spike theta phase for this part

            % work out the firing rate of the cluster in this part and add data to sdata
            sdata.celldata.(uci).(part_now).frate = numel(pspx) / part_duration;
            sdata.position.(part_now).duration = part_duration;         

% figure
% subplot(1,2,1)
% plot(ppox,ppoy,'k')
% hold on
% plot(pspx,pspy,'r.','MarkerSize',15)
% 
% subplot(1,2,2)
% plot(ppox,ppoy,'k')
% hold on
% plot(cspt,cspt,'r.','MarkerSize',15)
% 
% 
% keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Waveform data        
            if ~numel(pspx) % if there are no spikes there is no point in continuing
                continue
            end

            % get the waveforms for this cluster                
            waves2{1} = waves{1}(sindax,:);
            waves2{2} = waves{2}(sindax,:);
            waves2{3} = waves{3}(sindax,:);
            waves2{4} = waves{4}(sindax,:);
            for w = 1:4 % for every recording channel
                wav = double(waves2{w});
                ch = squeeze(mean(wav,1));
                chs = squeeze(std(wav,1));
                [maxval,maxt] = max(ch);
                [postminval,postmint] = min(ch(maxt:end));
                postmint = postmint + maxt - 1;
                width = postmint - maxt;
                width = width * (1000/50); 

                sdata.celldata.(uci).(part_now).waveform_mean{w} = ch;
                sdata.celldata.(uci).(part_now).waveform_stdv{w} = chs;
                sdata.celldata.(uci).(part_now).waveform_max(w) = maxval;
                sdata.celldata.(uci).(part_now).waveform_min(w) = postminval;
                sdata.celldata.(uci).(part_now).waveform_width(w) = width;
            end
            
            % calculate the signal to noise ratio of each channel - the root mean square of each signal / the root mean squared noise signal (cluster 0)
            rmsn = NaN(4,1);
            rmsn(1,1) = mean(rms(waves{1}(clu_identity == 0,:)));
            rmsn(2,1) = mean(rms(waves{2}(clu_identity == 0,:)));
            rmsn(3,1) = mean(rms(waves{3}(clu_identity == 0,:)));
            rmsn(4,1) = mean(rms(waves{4}(clu_identity == 0,:)));            
            rmss = NaN(4,1);
            rmss(1,1) = mean(rms(waves{1}(clu_indx,:)));
            rmss(2,1) = mean(rms(waves{2}(clu_indx,:)));
            rmss(3,1) = mean(rms(waves{3}(clu_indx,:)));
            rmss(4,1) = mean(rms(waves{4}(clu_indx,:)));                
            sdata.celldata.(uci).(part_now).channel_snr = rmss ./ rmsn; % accumulate data
            clear waves2 % clear all the waveform data, we don't need it again and it is quite large
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate ratemap and dwellmap    
            config.pos_ppm = pppm;
            config.spk_ppm = pspm;
            [ratemap,mapdata] = mapDATA4([ppox ppoy],[pspx pspy],config);
            dwellmap = mapdata.dwellmap;          
            
            % accumulate data
            sdata.celldata.(uci).(part_now).dwellmap = single(dwellmap); % add data to structure            
            sdata.celldata.(uci).(part_now).ratemap = single(ratemap); % add data to structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ratemap analyses          
            % ratemap analysis
            spatial_measures = getSPATinfo(ratemap,dwellmap);
            sdata.celldata.(uci).(part_now).spatial_measures = spatial_measures; % add data to structure          

            % place field analysis
            fieldd = getPFIELDS(ratemap,config.frcut,config.minfr,config.arcut);

            % accumulate data
            sdata.celldata.(uci).(part_now).field_count = uint8(length(fieldd.fields(:,1))); % add data to structure
            sdata.celldata.(uci).(part_now).field_data = fieldd; % add data to structure            

            % create autocorrelation
            automap = ndautoCORR(ratemap,ratemap,50);

            % autocorrelation analysis
            [grid_score,gdata] = gridSCORE(automap);

            % accumulate data
            sdata.celldata.(uci).(part_now).grid_autocorrelation = single(automap); % add data to structure
            sdata.celldata.(uci).(part_now).grid_score = grid_score; % add data to structure
            sdata.celldata.(uci).(part_now).grid_metrics = gdata;                        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overdispersion analysis     
            % convert position coordinates to map coordinates
            poxmap = mapdata.poxnew;
            poymap = mapdata.poynew;

            % generate smoothed spike histogram (instantaneous firing rate)
            [ppot2,uindx] = unique(ppot);
            poxmap2 = poxmap(uindx);
            poymap2 = poymap(uindx);
            pspt2 = unique(pspt);
            d2 = diff(ppot2)/2;

            edges = [ppot2(1)-d2(1); ppot2(1:end-1)+d2; ppot2(end)+d2(end)];
            edges(2:end) = edges(2:end)+eps(edges(2:end));
            shist = histcounts(pspt2,edges); % spike histogram
            fh = fspecial('gaussian',[1 config.over_smooth]); % boxcar filter
            shist = imfilter(shist,fh,'replicate','same'); % smooth histogram
            obsvals = shist(:) ./ 0.02; % convert to firing rate

            % calculate overdispersion
            expvals = ratemap(sub2ind(size(ratemap),round(poymap2),round(poxmap2)));
            zd = (obsvals - expvals) ./ sqrt(expvals);
            overd = nanstd(zd)^2;
            
            % accumulate data
            sdata.celldata.(uci).(part_now).over_dispersion = overd; % add data to structure            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HD analyses                 
            % Head direction analysis
            ai = linspace(0,2*pi,config.hd_bins)'; % angles for binning
            ai = ai(:);
            
            % HD density estimate
            hd_s = deg2rad(ppoh); % session head direction in radians
            hd_c = deg2rad(psph); % cell head direction in radians
            [hd1] = circ_ksdensity(hd_s,ai,[],config.hd_sigma); % the session head direction       
            [hd2] = circ_ksdensity(hd_c,ai,[],config.hd_sigma); % the cell's head direction
            hd1 = hd1 .* pos_tb; % convert session HD to time
            hd3 = hd2 ./ hd1; % calculate HD firing rate
            hd1 = hd1 ./ max(hd1); % normalise session hd
            hd3 = hd3 ./ max(hd3); % normalise cell hd
            hd3 = hd3(:);

            % head direction analyses
            rayleigh = circ_r(ai,hd3(:)); % rayleigh vector length
            mx2 = rad2deg(ai(hd3 == max(hd3))); % preferred angle (location of max frate)
            mn2 = rad2deg(circ_mean(ai,hd3)); % mean angle
            sd2 = rad2deg(circ_std(ai,hd3)); % std deviation angle
            
            % accumulate data
            sdata.celldata.(uci).(part_now).hd_density_session = hd1; % add data to structure
            sdata.celldata.(uci).(part_now).hd_density_cell = hd3; % add data to structure
            sdata.celldata.(uci).(part_now).hd_density_frate = hd_c; % add data to structure
            sdata.celldata.(uci).(part_now).hd_density_rayleigh = rayleigh; % add data to structure
            sdata.celldata.(uci).(part_now).hd_density_maximum = mx2; % add data to structure                
            sdata.celldata.(uci).(part_now).hd_density_mean = mn2; % add data to structure                
            sdata.celldata.(uci).(part_now).hd_density_stdev = sd2; % add data to structure                
                
            % HD histogram
            hd1 = hist(deg2rad(ppoh),config.hd_bins); % the session head direction   
            hd2 = hist(deg2rad(psph),config.hd_bins); % the cell's head direction
            fh = fspecial('average',[1 config.hd_boxcar]); % boxcar filter
            hd1 = imfilter(hd1,fh,'circular','same');            
            hd2 = imfilter(hd2,fh,'circular','same');            
            hd1 = hd1 .* pos_tb; % convert session HD to time
            hd3 = hd2 ./ hd1; % calculate HD firing rate
            hd1 = hd1 ./ max(hd1); % normalise session hd
            hd3 = hd3 ./ max(hd3); % normalise cell hd
            hd3 = hd3(:);
            
            % head direction analyses
            rayleigh = circ_r(ai,hd3); % rayleigh vector length
            mx2 = rad2deg(ai(hd3 == max(hd3))); % preferred angle (location of max frate)
            mn2 = rad2deg(circ_mean(ai,hd3)); % mean angle
            sd2 = rad2deg(circ_std(ai,hd3)); % std deviation angle
            
            % accumulate data
            sdata.celldata.(uci).(part_now).hd_session = hd1; % add data to structure
            sdata.celldata.(uci).(part_now).hd_cell = hd3; % add data to structure
            sdata.celldata.(uci).(part_now).hd_frate = hd_c; % add data to structure
            sdata.celldata.(uci).(part_now).hd_rayleigh = rayleigh; % add data to structure
            sdata.celldata.(uci).(part_now).hd_maximum = mx2; % add data to structure
            sdata.celldata.(uci).(part_now).hd_mean = mn2; % add data to structure                
            sdata.celldata.(uci).(part_now).hd_stdev = sd2; % add data to structure                
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike phase analyses              
            % Spikes vs time 
            bstvals = (min(ppot):config.time_bins:max(ppot)); % vector of time points at which we should calculate spike probability
            [bspikes,~] = histc(pspt,bstvals);
            [bsprobs,xibs] = ksdensity(pspt,bstvals);

            % accumulate data
            sdata.celldata.(uci).(part_now).spikes_time_histogram = [bstvals(:),bspikes(:)]; % add data to structure
            sdata.celldata.(uci).(part_now).spikes_time_ksdensity = [xibs(:),bsprobs(:)]; % add data to structure

            % bin the theta phase data
            ai = -pi:0.1:3*pi; 
            ai = ai(:);
            pspp2 = [pspp; (pspp+2*pi)];
            pspp2 = pspp2(:);
            yi = histc(pspp2,ai);
            phprobs = circ_ksdensity(pspp2,ai);
            
            % directional analyses on phase data
            rayleigh2r = circ_r(pspp(:)); % rayleigh vector length
            rayleigh2p = circ_rtest(pspp(:)); % rayleigh test for non-uniformity
            mx2p = ai(yi == max(yi)); % preferred angle (location of max frate)
            
            % example sine wave
            swav = cos(ai);
            swav = ((swav ./ max(swav))+abs(min((swav ./ max(swav))))).*max(yi);
            phmu = circ_mean(pspp);

            % accumulate data
            sdata.celldata.(uci).(part_now).phase_mean = phmu; % add data to structure
            sdata.celldata.(uci).(part_now).phase_r = rayleigh2r; % add data to structure
            sdata.celldata.(uci).(part_now).phase_p = rayleigh2p; % add data to structure            
            sdata.celldata.(uci).(part_now).phase_max = mx2p(1); % add data to structure
            sdata.celldata.(uci).(part_now).spike_phase_ideal = swav(:); % add data to structure
            sdata.celldata.(uci).(part_now).spike_phase_binned = yi(:); % add data to structure            
            sdata.celldata.(uci).(part_now).spike_phase_ksdensity = phprobs(:); % add data to structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Theta and bursting analyses   
% Theta index from Mizuseki et al. (2012) Activity Dynamics and Behavioral Correlates of CA3 and CA1 Hippocampal Pyramidal Neurons 
% Burst index also from Mizuseki et al. (2012) Activity Dynamics and Behavioral Correlates of CA3 and CA1 Hippocampal Pyramidal Neurons 
% Theta skipping index from Brandon et al. (2013) Segregation of cortical head direction cell assemblies on alternating theta cycles - theta_skip is bound between -1 and 1, and higher values indicate more theta cycle skipping.
            % burst index
            isis = diff(pspt)*1000; % interspike intervals in ms
            bindx = zeros(size(pspt));
            bindx([isis; NaN] < 6 | [NaN; isis] < 6) = 1; % bindx is an index of all spikes sharing an isi less than 6ms
            bindx = logical(bindx);
            burst_index = sum(bindx) / numel(pspt);
            
            sts = regionprops(bindx,'Area');
            burst_length_median = nanmedian([sts.Area].');
            burst_length_mean = nanmean([sts.Area].');

            % accumulate data
            sdata.celldata.(uci).(part_now).burst_index = burst_index;
            sdata.celldata.(uci).(part_now).burst_length_median = burst_length_median;
            sdata.celldata.(uci).(part_now).burst_length_mean = burst_length_mean;

            % Spike autocorrelation - theta analaysis           
            % 250ms autocorrelation
            [tms1,corrdata1] = spikeINTERVALS_v2(pspt,250,1);  
            sdata.celldata.(uci).(part_now).theta_train = [corrdata1(:) tms1(:)]; % accumulate data
            
            % 400ms autocorrelation    
            [tms3,corrdata3] = spikeINTERVALS_v2(pspt,400,1);  
            corrdata3(tms3 == 0) = 0;         
            uni_corr = corrdata3 ./ nanmax(corrdata3(:)); % we want the sum of probability to be unity
            
            [b,a] = butter(4,[5 12]/(1000/2)); % Generate 4th order butterworth filter coefficients
            uni_corr_filt = filtfilt(b,a,uni_corr); % Apply filter to data using zero-phase filtering  
            hilbert_corr = hilbert(uni_corr_filt); % calculate hilbert transform of the filtered signal    
            theta_index = nanmean(abs(hilbert_corr)); % take the mean of the abs of the hilbert transform (the envelope of the signal or the instantaneous amplitude)
            theta_index = theta_index * 1000; % scale the index up to something we actually care about
            rindx = knnsearch(tms3(:),[-65 65]');
            theta_ratio = nansum([corrdata3(rindx(1)) corrdata3(rindx(2))]) ./ nansum([corrdata3(1) corrdata3(end)]); % theta ratio: the ratio of events at the points [-65 and 65] ms to events at the points [-400 and 400] ms

            % theta skipping ratio
            [tms_skip,prob_skip,~] = spikeINTERVALS_v2(pspt,400,10);    
            tms_skip = tms_skip ./ 1000;
            [xData,yData] = prepareCurveData(tms_skip,prob_skip);
            ft = fittype( '(a1*(cos(w*x)+1) + a2*(cos(0.5*w*x)+1) + b) * (exp(-abs(x)/t1)) + (c.*exp(-x^2/(t2^2)))', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions('Method','NonlinearLeastSquares');
            opts.Algorithm = 'Trust-Region';
            opts.Display = 'Off';
            m = nanmax(prob_skip);
            opts.Lower = [0 0 0 -m 0 0 pi*10];
            opts.StartPoint = [0.865 0.042 0.931 0.674 0.567 0.226 0.988];
            opts.Upper = [m m m m 5 0.05 pi*18];
            [fitresult,~] = fit(xData,yData,ft,opts);
            p1x = (2*pi)/fitresult.w;
            p1 = fitresult(p1x);
            p2x = (4*pi)/fitresult.w;
            p2 = fitresult(p2x);
            theta_skip = (p2-p1) ./ max([p2 p1]);                 
            
            % accumulate data        
            sdata.celldata.(uci).(part_now).theta_train_long = [corrdata3(:) tms3(:) uni_corr_filt(:)];
            sdata.celldata.(uci).(part_now).theta_index = theta_index;
            sdata.celldata.(uci).(part_now).theta_ratio = theta_ratio;
            sdata.celldata.(uci).(part_now).theta_skip = theta_skip;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speed modulation analyses      
% from Kropff, Carmichael, Moser and Moser (2015) Speed cells in the medial entorhinal cortex
            % paper says they used a gaussian filter over 250ms, but samples are every 20ms... so I have used 13 samples or 260ms instead
            fh = fspecial('gaussian',[1 13]); % boxcar filter
            shist = imfilter(shist,fh,'replicate','same'); % smooth histogram
            shist = shist ./ 0.02; % convert to firing rate

            % for figures we want to bin firing rate relative to speed, in 2cm/s increments, up to 50cm/s, excluding 0-2cm/s
            bedges = 0.02:0.02:0.5;
            bvals = zeros(numel(bedges)-1,1);
            for bb = 1:numel(bvals)
                bindx = ppov(uindx) > bedges(bb) & ppov(uindx) <= bedges(bb+1);
                mval = nanmean(shist(bindx));
                bvals(bb) = mval;
            end
            xvals = mean([bedges NaN; NaN bedges],1);
            xvals = xvals(2:end-1);
            xvals = xvals(:)*100;
            bvals = bvals(:);

            % statistics
            sscore = corr(xvals,bvals,'type','Pearson','rows','pairwise'); % speed score is the correlation between instantaneous firing rate and speed
            P = polyfit(xvals,bvals,1);
            bvals2 = (bvals - P(2)) / (P(1) * 0.50);
            
            % accumulate data
            sdata.celldata.(uci).(part_now).speed_score = sscore;
            sdata.celldata.(uci).(part_now).speed_slope = P(1);
            sdata.celldata.(uci).(part_now).speed_intercept = P(2);
            sdata.celldata.(uci).(part_now).speed_curve = [xvals bvals];
            sdata.celldata.(uci).(part_now).speed_curve_normalised = [xvals bvals2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Refractory period analyses        
% Refractory period contamination from Navratilova and McNaughton (2016) Grids from bands, or bands from grids? An examination of the effects of single unit contamination on grid cell firing fields
% and Fee, Mitra, Kleinfeld (1996) Automatic sorting of multiple unit neuronal signals in the presence of anisotropic and non-Gaussian variability
            % Spike autocorrelation - refractory period analysis
            [tms2,corrdata2] = spikeINTERVALS_v2(pspt,25,0.5);        

            % Calculate refractory contamination
            half_spike = corrdata2(tms2 >= 0); % take only the positive side of spike autocorrelogram
            half_time = tms2(tms2 >= 0); % take only the positive side of spike autocorrelogram
            tau_tot = tau_r - tau_c;
            Ns = numel(ppox); % number of cluster spikes
            lambda = Ns/part_duration;  % mean firing rate for cluster 
            RPV = sum((half_spike ~= 0) & (half_time <= tau_r)); % get the number of refractory period violations

            % get Poisson confidence interval on number of expected RPVs
            conf_int = 95; % percent confidence interval
            [~,interval] = poissfit(RPV,(100-conf_int)/100); 

            % convert contamination from number of RPVs to a percentage of spikes
            RPVs = [RPV interval(1) interval(2)];
            cont_bounds = NaN(1,3);
            for i = 1:length(RPVs)
                RPVnow = RPVs(i);
                RPVT = 2 * tau_tot * Ns; % total amount of time in which an RPV could occur = the usable refractory period (tau_tot) around each spike (*Ns) and on each side of the spike (*2)
                RPV_lambda = RPVnow / RPVT; % rate of RPV occurence
                p =  RPV_lambda / lambda; % estimate of % contamination of cluster

                % force p to be a real number in [0 1]
                if isnan(p)
                    p = 0; 
                end
                if p > 1
                    p = 1; 
                end   % if p > 1
                cont_bounds(i) = p;
            end
            censored_estimate = (tau_c/1000) * spik_count / part_duration; % estimate the fraction of spikes not detected because of the system lockout

            % accumulate data
            sdata.settings.tau_r = tau_r;
            sdata.settings.tau_c = tau_c;
            sdata.celldata.(uci).(part_now).refractory_period = [corrdata2(:) tms2(:)];
            sdata.celldata.(uci).(part_now).refractory_violations = RPV; % add data to structure
            sdata.celldata.(uci).(part_now).refractory_contamination = cont_bounds(1); % add data to structure
            sdata.celldata.(uci).(part_now).refractory_contamination_95_lower_upper_bounds = cont_bounds(2:3); % add data to structure
            sdata.celldata.(uci).(part_now).refractory_censoring = cont_bounds(2:3); % add data to structure                
            sdata.celldata.(uci).(part_now).censored_estimate = censored_estimate; % add data to structure                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell type/identity
            [ctype,ctypen,ctypeb] = getCELLTYPE(sdata,uci,part_now);
            sdata.celldata.(uci).(part_now).cell_type = ctype; % add data to structure                
            sdata.celldata.(uci).(part_now).cell_type_num = ctypen; % add data to structure   
            sdata.celldata.(uci).(part_now).cell_type_bin = ctypeb; % add data to structure               
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures
            %% Overall summary figure
            if run_figCLUST
                figCLUST(sdata,uci,part_now,fig_vis,save_fig); % overall figure with ratemaps etc
            end 
        end 
        
        %% Parts summary figure
        if run_figPARTS
            if nparts > 1 % if there is more than one part (otherwise this figure is not useful)
                figPARTS(sdata,uci,fig_vis,save_fig); % figure with ratemaps for each part
            end 
        end
        disp(sprintf('\b %.f%%',cc/length(clus)*100));
    end

    %% Cluster cross-correlation figure
    if run_figCROSS
        if sum(clus ~= 0) > 1
            [~,~,~] = mkdir([pwd '\klustest\' sdata.settings.combined_name '\figures\']);
            figfile = [pwd '\klustest\' sdata.settings.combined_name '\figures\' sdata.combined_name '_E' num2str(tet) '_cross_correlogram.png'];            
            figCROSS(sdata,tet,figfile,fig_vis); % overall figure with ratemaps etc
        end
    end

    %% Cluster space figure
    if run_figCSPACE
        [~,~,~] = mkdir([pwd '\klustest\' sdata.settings.combined_name '\figures\']);
        figfile = [pwd '\klustest\' sdata.settings.combined_name '\figures\' sdata.combined_name '_E' num2str(tet) '_cluster_space.png'];
        figCSPACE(sdata.elecdata.(tetstr).fetdata,figfile,fig_vis); % overall figure with ratemaps etc
    end
end       

%% Save the session data structure
save([pwd '\klustest\' cname '\' cname '_sdata.mat'],'sdata','-v7.3'); % save session data
save([pwd '\klustest\' cname '\' cname '_mtint.mat'],'mtint','-v7.3'); % save mtint file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finish up
toc1 = toc/60;
disp(sprintf('klustest has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',[pwd '\klustest\' cname '\figures'],' &'');">','figures folder','</a>'])
disp('-------------------------------------------------------------------------------------');

% send an email if requested
if alert_me
    sub = sprintf('klustest finished at: %s',datestr(now));
    [~,compname] = system('hostname');
    msg1 = sprintf('%s started running klustest at: %s \nIt took %.f seconds or %.2f minutes',compname,d_start,toc,toc/60);
    MailMe(alert_email,sub,msg1)
end





























