function klustest(tetrodes,clusters,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%klustest  cluster analysis function
% This function utilises many minor functions to both analyse and plot cluster data generated 
% in Tint. If requested, it will output a figure for each cluster, the cluster space of each
% tetrode, the cross-correlations of every cluster on a tetrode and a session data structure (sdata.mat)
% It will also generate an mtint file (mtint.mat) containing all the tetrode and cluster info.
% klustest(tetrodes,clusters,rname,cname)
%
% USAGE:
%         klustest(tetrodes,clusters,rname,cname)
%
% INPUT:
%         tetrodes - (default = 1:16) the tetrodes to run on in a vector format (i.e. [1 2 3 4 5 6] would run on tetrodes 1 to 6)
%         clusters - (default = 0) the clusters to run on, set this to 0 to run on all clusters
%         rname - the rat name/number as a string - this will be used in the sdata structure so its important to give this
%         cname - (optional) function will automatically look for klustakwiked files named 'kwiktint', but you can change this here (i.e if you chose a custom output name in kwiktint)
%
% EXAMPLES:
%
% See also: KWIKTINT

% HISTORY:
% version 1.0.0, Release 05/08/16 created from an older version
% version 1.1.0, Release 08/08/16 modified readDACQDATA to output waveforms and added post processing of this info
% version 1.2.0, Release 08/08/16 readDACQDATA working, added postprocessing of mtint
% version 1.3.0, Release 09/08/16 main plots working
% version 1.4.0, Release 10/08/16 options added for different styles of plot
% version 2.0.0, Release 11/08/16 added theta autocorrelation sine estimate (O'Mara 2015)
% version 2.1.0, Release 12/08/16 added refractory period violations (Navratilova and McNaughton (2016)
% version 2.1.1, Release 13/08/16 added saving data, figures
% version 2.2.0, Release 14/08/16 added cell type identification
% version 3.0.0, Release 15/08/16 cluster quality assessment added in clusterQUALITY
% version 3.1.0, Release 16/08/16 added cluster space figures
% version 3.1.1, Release 17/08/16 fixed phase plot and phase map
% version 3.2.0, Release 18/08/16 modified clusterQUALITY to deal better with missing channels
% version 3.3.0, Release 18/08/16 added the option to ignore position data, fixed bug with waveform plot
% version 3.3.1, Release 19/08/16 make it so none of the position figures are made if there is no position data, this should be faster
% version 3.4.0, Release 22/08/16 fixed issues with cluster space plot, concentrate on first feature, added legend
% version 3.5.0, Release 23/08/16 added cluster space subplot to cell figure, for Ele, this uses the first clustering feature, the 1st and 2nd highest amplitude channel
% version 4.0.0, Release 25/06/16 added tetrode input to readalldacqdata, this prevents it trying to open unnecessary files, fixed bug in the detection of dead/lfp channels
% version 4.1.0, Release 25/08/16 fixed a minor bug in histograms when there is only one spike
% version 5.0.0, Release 22/11/16 adapted from KlusterAnalysisTINT
% version 5.1.0, Release 30/03/17 started editing the newer version with more emphasis on sdata structure
% version 6.0.0, Release 30/03/17 adapted to new partitioning system using a part_config cell array, function can now split based on any inputs
% version 7.0.0, Release 31/03/17 started moving figure work into subfunctions to contain them, each uses the sdata structure
% version 7.1.0, Release 01/04/17 finished most of the figure subfunctions
% version 8.0.0, Release 02/04/17 added kluspartfig, which generates a summary figure containing all the parts, and correlations between them
% version 9.0.0, Release 05/04/17 fixed multiple problems in clusterQUALITY, including wrong lratios and reliance on kwiktint files, now using clusterQUALITY_v2
% version 10.0.0, Release 13/04/17 added ability to extract parent filenames from .cut file, meaning the function only needs the kwiktint out name
% version 10.1.0, Release 19/04/17 fixed head direction plot and analyses
% version 10.2.0, Release 20/04/17 changed from using mapDATA to mapDATA_v2, latter allows for maps with same limits and converted position data necessary for overdispersion
% version 10.3.0, Release 20/04/17 added overdispersion calculation (my own)
% version 10.4.0, Release 05/05/17 added comments, added easier naming of part_config fields
% version 10.5.0, Release 10/05/17 fixed bugs with interval times in part_config, added figKEYS to handle interval times, added exceptions to figCLUST and overdispersion to handle overlapping trials
% version 10.6.0, Release 12/05/17 fixed bug in getDACQDATA where keypress times were not accumulated in time, fixed bug in partSESS where incorrect interval vector was being loaded
% version 11.0.0, Release 12/05/17 added ability to send emails when completed
% version 11.1.0, Release 17/05/17 added possibility to calculate shuffled spatial measures
% version 11.2.0, Release 19/05/17 replaced read_key with saveKEY, this automatically saves a text file version of the .inp files as it loads them
% version 12.0.0, Release 19/05/17 replaced getPARTconfig with getPARTconfig_v2, this saves a much nicer text file which is easier to edit
% version 12.1.0, Release 07/06/17 fixed a problem in getDACQDATA where it would try to run on all tetrodes even if a subset is specified
% version 12.2.0, Release 22/06/17 fixed error in time allocation
% version 12.3.0, Release 10/07/17 changed part_config specification so that it is done in a loop, added information required by contest and klustest3
% version 13.0.0, Release 25/07/17 added variable pixel ratio capability (big job!) getDACQDATA now uses getDACQDATAHEADERS instead of the horible key_value functions
% version 13.1.0, Release 01/11/17 fixed bug where LFP and session duration did not match
% version 14.0.0, Release 01/11/17 replaced celltype with getCELLTYPE, replaced mapDATA_v3 with mapDATA4, replaced GridAnalaysis with gridSCORE2, replaced AutoCorr with ndautoCORR
% version 15.0.0, Release 04/08/18 overhauled, code cleaned up, data maintained in cm throughout
% version 15.1.0, Release 10/08/18 recoded overdispersion and speed score analysis
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% I can extract the rat name from the file path because I know my folder structure, if this is not the case for you, comment out the next line
pname = pwd; sindx=strfind(pname,'\'); rname=pname(sindx(end-1)+1:sindx(end-1)+3);

% deal with the other variables
inps = {'tetrodes','clusters','cname','rname'};
vals = {'1:16','0','''kwiktint''','''000'''};
for ff = 1:length(inps)
    if ~exist(inps{ff},'var')
        eval([inps{ff} '=' vals{ff} ';']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIAL SETTINGS / INPUTS
    % Specify partition configuration settings
    part_names = {'square1','lattice','square2'}; % the names you want the outputs to be saved as, these cannot start with a number: i.e. 'session1'
    % method of partition, corresponding to each of the names above: 1 = combine everything, 2 = take a recording session, 3 = use digital inputs
    part_methods = [2 2 2 2]; % i.e. if part_methods=[2 2 2 1], we will have 4 parts, the first 3 correspond to some combination of recording sessions (there can be multiple ones) and the last one will include all data
    % cell array of vectors indicating which intervals to include in each partition: if method = 1 this does nothing, if method = 2 this should specify which recording sessions to include, if method = 3 this should specify which digital input pairs to include (inf = use all)
    part_intervals = {1 2 3 4}; % i.e. if part_methods=[2 2 2], then if part_intervals={1 2 3}, rec 1 will go in part 1, rec 2 in part 2 and rec 3 in part 3    OR     if part_methods=[2 2 2 2 2], then if part_intervals={1 [2 4 5] 3}, rec 1 will go in part 1, rec 2,4 and 5 in part 2 and rec 3 in part 3
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
    part_config.FDATA.part_names = part_names; % the manually defined part_names
    part_config.FDATA.out_name = cname; % the name to use when finding kwiktint files and saving new files, default is kwiktint
    part_config.FDATA.rat_name = rname; % the rat number is a useful thing to include here as it cannot be retrieved anywhere else
    part_config.FDATA.recording_times = []; % starts empty, will contain the start and end times for every recording session
    % part_config.FDATA.enames = {'circ','cbox','ecbox'}; % for contest: these environment names correspond to etype
    % part_config.FDATA.pratio = [395 446 446 446 446 446 446 446 446 446 395]; % recording session pixel ratios, if this is not empty these will be used instead of the dacq data file values. This may not correspond to your parts unless you are partitioning by method 2

%     % Specify partition configuration settings
%     part_names = {'square1','lattice','square2'}; % the names you want the outputs to be saved as, these cannot start with a number: i.e. 'session1'
%     % method of partition, corresponding to each of the names above: 1 = combine everything, 2 = take a recording session, 3 = use digital inputs
%     part_methods = [2 2 2]; % i.e. if part_methods=[2 2 2 1], we will have 4 parts, the first 3 correspond to some combination of recording sessions (there can be multiple ones) and the last one will include all data
%     % cell array of vectors indicating which intervals to include in each partition: if method = 1 this does nothing, if method = 2 this should specify which recording sessions to include, if method = 3 this should specify which digital input pairs to include (inf = use all)
%     part_intervals = {1 2 3}; % i.e. if part_methods=[2 2 2], then if part_intervals={1 2 3}, rec 1 will go in part 1, rec 2 in part 2 and rec 3 in part 3    OR     if part_methods=[2 2 2 2 2], then if part_intervals={1 [2 4 5] 3}, rec 1 will go in part 1, rec 2,4 and 5 in part 2 and rec 3 in part 3
%     dimensions = [2 3 2]; % for klustest3
%     % part_etypes = [1 3 2 2 2 2 2 2 2 3 1]; % for contest: environment order: 0 = ignore, 1 = circle, 2 = contextbox, 3 = empty contextbox
%     % part_rtypes = [0 0 0 90 0 270 0 270 180 0 0]; % for contest: clockwise rotation angle
%     for pp = 1:length(part_names)
%         part_config.(part_names{pp}).method = part_methods(pp);
%         part_config.(part_names{pp}).intervals = cell2mat(part_intervals(pp));
%         part_config.(part_names{pp}).inputs = []; % starts empty, this will contain keypress-time pairs, used if method = 3
%         part_config.(part_names{pp}).times = []; % starts empty, this will specify the start and end time(s) in the combined file for this partition
%         part_config.(part_names{pp}).dimensions = dimensions(pp); % for klustest3
%     %     part_config.(part_names{pp}).etype = part_etypes(pp); % for contest
%     %     part_config.(part_names{pp}).rtype = part_rtypes(pp); % for contest
%     end
%     % any fields under FDATA are saved along with part_config, but are not used for partitioning
%     part_config.FDATA.part_names = part_names; % the manually defined part_names
%     part_config.FDATA.out_name = cname; % the name to use when finding kwiktint files and saving new files, default is kwiktint
%     part_config.FDATA.rat_name = rname; % the rat number is a useful thing to include here as it cannot be retrieved anywhere else
%     part_config.FDATA.recording_times = []; % starts empty, will contain the start and end times for every recording session
%     % part_config.FDATA.enames = {'circ','cbox','ecbox'}; % for contest: these environment names correspond to etype
%     % part_config.FDATA.pratio = [395 446 446 446 446 446 446 446 446 446 395]; % recording session pixel ratios, if this is not empty these will be used instead of the dacq data file values. This may not correspond to your parts unless you are partitioning by method 2

% overrides
pconfig_override = 1; % set to 1 if you want to ignore and overwrite an existing part_config, set to 2 to run with current part_config settings without overwriting anything
maintain_mtint = 0; % set to 1 to save/load mtint in the base workspace, this saves time when running the function mutliple times (for instance in debugging) but should otherwise be set to 0
mtint_override = 0; % set this to 1 to force the production of a new mtint file (do this if you reclustered in tint for instance)
lattice_pos = 0; % set to 1 to load 3D trajectory produced by lattice instead of dacqUSB trajectory

% map settings
config.rmethod = 'gaussian'; % (default 'nearest') the mapping approach to use, either 'nearest','gaussian','adaptive','KDE'
config.map_padd = 2; % (default 2) the number of bins to pad spatial maps with
config.bin_size = 2; % (default 2) bin size in cm for calculating the rate map (make sure data is in cm or pm/sm values are given)
config.map_sigma = 1.5; % (default 1.5) used by nearest and KDE method, sigma of gaussian to use when smoothing traditional dwell and spike maps, or used as bandwidth of KDE
config.min_dwell = 0.01; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
config.g_sigma = 5; % (default 10) only used by gaussian method - sigma of gaussian used to weight position point and spike point distances, bigger means smoother maps
config.min_dist = 5; % (default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
config.max_dist = 20; % (default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
config.srate = 50; % (default 50) sampling rate of data in Hz, used to calculate time

% field settings
config.frcut = 0.5; % relative minimum firing rate (% of ratemap max) to be considered a field
config.arcut = 9; % minimum number of pixels to be considered a field, typically people use 9 contiguous pixels
config.minfr = 1; % (Hz) absolute minimum cutoff firing rate to be considered a field

% spike plot settings
config.over_smooth = 13; % number of position data points over which to smooth instantaneous firing rate when calculating overdispersion
config.time_bins = 2; % (default 2s) time window over which to compute the spike vs time plot
            
% HD settings
config.hd_type = 'histogram'; % (default 'histogram') enter 'density' for a circular kolmogorov smirnov density estimate plot, enter 'histogram' for the traditional histogram polar plot
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

% optimization settings
config.waveform_analyses = 1;
config.cluster_analyses = 1;
config.temporal_analyses = 1;
config.field_analyses = 1;
config.grid_analyses = 1;
config.hd_analyses = 1;
config.overdispersion_analyses = 1;
config.spike_isi_analyses = 1;
config.spike_autocorrelation_analyses = 1;
config.spike_phase_analyses = 1;
config.speed_analyses = 1;
config.cell_analyses = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start analysis
% get the function's name (people may rename it from klustest)
stk = dbstack;
function_name = stk.name;

% starting messages
tic;
disp('----------------------------------------------------------------------------');
disp(sprintf('Running %s at %s...',function_name,datestr(now,'yyyy-mm-dd-HH-MM-SS')))

% prepare sdata
sdata = struct; % create an empty structure - together with the mtint file this will hold all of the session data

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
disp(sprintf('\t...rat name: %s',rname))    

% collect the data
sdata.combined_name = cname; % add data to structure
sdata.rat_num = rname;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get session data
% get session names from cut file
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

% check which tetrodes actually exist
disp(sprintf('Assessing data...'))
[tetrodes,mvalue] = getTRODES(snames,tetrodes);
if isempty(mvalue)
    disp(sprintf('\t...tetrodes: %s accounted for',mat2str(tetrodes)))
else
    disp(sprintf('\t...WARNING: tetrodes: %s are incomplete and will be skipped',mat2str(mvalue)))
end
disp(sprintf('\t...done'))

% collect data
sdata.session_names = snames; % add data to structure
sdata.num_sessions = nsess;
sdata.config = config;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read all dacq data or load it
disp(sprintf('Fetching DACQ data...'));
if ~maintain_mtint || mtint_override
    clear mtint
end

mname = ['klustest\' cname '\' cname '_mtint.mat'];
if any(strcmp(evalin('base','who'),'mtint')) && ~mtint_override && maintain_mtint
    disp(sprintf('\t...using mtint held in memory'));
    mtint = evalin('base','mtint');
elseif ~exist(mname,'file') || mtint_override
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
    save(mname,'-struct','mtint','-v7.3');
    
    if maintain_mtint
        assignin('base','mtint',mtint); % leave mtint in base workspace
    end    

    sdata.date = mtint.header_date;
    part_config.FDATA.recording_times = mtint.pos.trial_duration;
    pos = mtint.pos;
    tetrode = mtint.tetrode;
    lfp = mtint.lfp;
else
    disp(sprintf('\t...using saved mtint'));
    
    % I have changed the mtint to save as variables rather than one structure (upon loading it is always a structure though)
    % to avoid conflicts with earlier versions of klustest we have to check here what format the mtint is in
    try
        warning('off','MATLAB:load:variableNotFound')
        load(mname,'header_date');
        sdata.date = header_date;
        load(mname,'pos');
        part_config.FDATA.recording_times = pos.trial_duration;
        load(mname,'tetrode');      
        load(mname,'lfp');      
    catch
        load(mname);
        sdata.date = mtint.header_date;
        part_config.FDATA.recording_times = mtint.pos.trial_duration;
        pos = mtint.pos;
        tetrode = mtint.tetrode;
        lfp = mtint.lfp;
    end
    warning('on','all'); 
end
disp(sprintf('\t...done'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% partitioning settings
disp(sprintf('Preparing partitions...'))
% determine partition time values
disp(sprintf('\t...finding part_times'));
part_config = partSESS(part_config,mname); % process the part_config to get the start and end time of each partition
getPARTconfig_v2(pconfig_name,part_config); % save the current part_config cell array
nparts = length(fieldnames(part_config))-1;
disp(sprintf('\t...done'));
sdata.part_config = part_config;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get some variables that will be important later
disp(sprintf('Extracting initial data...'))
disp(sprintf('\t...recording sessions: %d',numel(snames)))

% get the total length of the session
duration = pos.total_duration;
disp(sprintf('\t...total session time: %ds',duration))
sdata.session_duration = duration; % add data to structure

% get the position data (dacqUSB data) for the whole session
if lattice_pos
    disp(sprintf('\t...loading reconstructed trajectory instead of dacqUSB trajectory'))
    [pox,poy,~] = getLATTICEpath(snames,pos);
    pot = double(pos.ts); % extract the time stamp of each position value
    pov = double(pos.speed(:)); % the running speed throughout the session
    poh = double(pos.dir(:,1)); % HD data
    ppm = ones(size(pot)).*1000; % pixels per metre     
else
    position = pos.xy_pixels;	
    posx = position(:,1); % extract just the x coordinates
    posy = position(:,2); % extract just the y coordinates
    pox = double(posx);
    poy = double(posy);
    pot = double(pos.ts); % extract the time stamp of each position value
    pov = double(pos.speed(:)); % the running speed throughout the session
    poh = double(pos.dir(:,1)); % HD data
    ppm = {pos.header.pixels_per_metre_vec}.'; ppm = double(vertcat(ppm{:,:})); % pixels per metre
end

% convert the data to cm
pox = pox ./ ppm .* 100; % convert position data to cm
poy = poy ./ ppm .* 100; % convert position data to cm
    
% adjust the data
com_min_x = min(pox);
com_min_y = min(-poy);
pox = pox - com_min_x;
poy = -poy - com_min_y;

% accumulate data
sdata.pox = pox;
sdata.poy = poy;
sdata.pot = pot; 
sdata.pov = pov;
sdata.ppm = ppm;
poss1 = numel(pox);
leds = size(pos.led_pos,2);

% get the position data sampling rate (should be 50hz) or 0.05s
samp_rate_hz = pos.header(1).sample_rate_num;
samp_rate_hz = samp_rate_hz(1,1);
pos_tb = 1 / samp_rate_hz;
config.srate = samp_rate_hz;

% display results
disp(sprintf('\t...positions read: %d',poss1));
disp(sprintf('\t...tracking LEDs: %d',leds));
disp(sprintf('\t...median pixel ratio: %dppm',nanmedian(ppm)));
disp(sprintf('\t...sample rate: %dHz (%.2fs)',samp_rate_hz,pos_tb));

% accumulate data
sdata.position_timebase = pos_tb;
sdata.position_srate_hz = samp_rate_hz;
sdata.leds = leds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve and analyse data
disp(sprintf('Analysing lfp data...'))
Fs = lfp(1).Fs(1,1);
lfp = double(lfp(1).lfp(:,1));
lfptime = (0:length(lfp)-1)'/Fs; % make a vector for time
disp(sprintf('\t...samples read: %d',numel(lfp)));
disp(sprintf('\t...sample rate: %dHz',Fs));
if round(max(lfptime)) ~= round(duration)
    disp(sprintf('\tWARNING: lfp time %.2f does not match session duration %.2f...',max(lfptime),duration));
end

% filter lfp to get theta
cutOffFreq = [4 12]; % lower and upper cutoff frequencies in Hz
[b,a] = butter(4,cutOffFreq/(Fs/2)); % Generate 4th order butterworth filter coefficients
lfpfilt = filtfilt(b,a,lfp); % Apply filter to data using zero-phase filtering
lfp_theta = lfpfilt;
lfp_t = lfptime;
disp(sprintf('\t...done'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve and analyse data
disp(sprintf('Analysing tetrode data...'))

%% For every available tetrode
for e = 1:length(tetrodes) 
    tet = tetrodes(e); % tet = the current tetrode
    tetstr = ['t' num2str(tet)];
    disp(sprintf('\tLoading tetrode %d...',tet));
    
    pos_samps = ceil(tetrode(tet).ts * 50);    
    spik_count = tetrode(tet).nspike_cut; % retrieve the number of spikes recorded on this tetrode
    sdata.spike_count{tet} = spik_count;
    
    % get a vector of clusters we want to analyse
    if clusters == 0
        clus = unique(tetrode(tet).cut); % find the clusters logged for this tetrode
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
    disp(sprintf('\t\t\t\t...starting analysis'));
    
    % calculate the quality of clusters in this session    
    if config.cluster_analyses
        disp(sprintf('\t\t\t\t...getting cluster quality'));    
        fetname = [sdata.combined_name '.fet.' num2str(tet)];
        cutname = [sdata.combined_name '_'  num2str(tet) '.cut'];
        fetdata = clusterQUALITY_v2(fetname,cutname);        
        sdata.(tetstr).fetdata = fetdata;
        sdata.(tetstr).clusters = clus;
        clear fetdata
    end
    
    % get the channel waveforms for this tetrode    
    if config.waveform_analyses
        disp(sprintf('\t\t\t\t...getting waveforms'));        
        waves = cell(1,4);
        for ggnow = 1:length(sdata.session_names)
            fnamen = sdata.session_names{ggnow};
            [~,c1,c2,c3,c4] = getspikes([fnamen,'.',num2str(tet)]);   
            waves{1} = [waves{1}; c1];
            waves{2} = [waves{2}; c2];
            waves{3} = [waves{3}; c3];
            waves{4} = [waves{4}; c4];  
        end
    end
    
    % get the spike times for this tetrode
    spiketime = tetrode(tet).ts;
     
%% For every detected cluster    
    dispPROGRESS(0); % initialise percentage progress

    for cc = 1:length(clus) % for every cluster
        clu = clus(cc); % clu = the current cluster
        uci = ['r' rname '_' sdata.date '_t' num2str(tet) '_c' num2str(clu)]; % unique cell identifier - a string with [rat name, session date, electrode, cluster];
                      
        % get some vectors that we can use to sort data
        clu_identity = tetrode(tet).cut; % clu_assign is a vector of numbers, one for each spike, each number corresponds to a cluster
        clu_indx = (clu_identity == clu); % logical vector, one number per spike, 1 if in cluster, 0 if not
        n_spikes = length(find(clu_identity == clu));
        pos_identity = pos_samps; % pos_identity is a vector of numbers, one for each spike, each number corresponds to a position data point
            
        % sort out spike data
        sdata.(uci).spike_times = spiketime(clu_indx); % the time point for every spike in this cluster
        
        % if this is the noise cluster don't continue any further
        if ~clu      
            dispPROGRESS(cc,1:length(clus)); % display progress percentage            
            continue
        end
        
%% For each part  
        part_names = fieldnames(part_config);
        for pp = 1:nparts % for every partition
            part_now = part_names{pp}; % the name of the current part
            part_times = part_config.(part_now).times; % the time pairs (intervals) corresponding to this part

            % find what data falls into the intervals associated with this part
            sindax = zeros(size(spiketime));
            pindax = zeros(size(pox));
            for ii = 1:length(part_times(:,1)) % for every interval
                sindax = logical(sindax + (spiketime > part_times(ii,1) & spiketime < part_times(ii,2) & clu_indx)); % logical array, same length as spiketime, 1 when spike is in this cluster and is in this part
                pindax = logical(pindax + (pot > part_times(ii,1) & pot < part_times(ii,2))); % logical array, same length as pox, 1 means position sample falls into this part
            end 
            part_duration = sum(pindax(:))*pos_tb;            
            
            % extract spike and position data
            ppox = pox(pindax); % pos x
            ppoy = poy(pindax); % pos y
            ppot = pot(pindax); % pos time
            ppov = pov(pindax); % pos running speed
            ppoh = poh(pindax); % pos HD
            pppm = ppm(pindax); % pixel ratio
            
            pspt = spiketime(sindax);
            sidx = knnsearch(pot,pspt);
            pspx = pox(sidx); % spike x
            pspy = poy(sidx); % spike y
            psph = poh(sidx); % direction
            pspv = pov(sidx); % running speed
            pspm = ppm(sidx); % pixel ratio

            % assign data to arrays
            sdata.(part_now).pox = single(ppox);
            sdata.(part_now).poy = single(ppoy);
            sdata.(part_now).pot = single(ppot); 
            sdata.(part_now).poh = single(ppoh); 
            sdata.(part_now).pov = single(ppov);
            sdata.(part_now).ppm = single(pppm);            
            sdata.(uci).(part_now).spx = single(pspx);
            sdata.(uci).(part_now).spy = single(pspy);
            sdata.(uci).(part_now).spt = single(pspt);
            sdata.(uci).(part_now).spv = single(pspv);
            sdata.(uci).(part_now).sph = single(psph);     
            sdata.(uci).(part_now).spm = single(pspm);     
            
            % work out the firing rate of the cluster in this part
            pfrate = numel(pspx) / part_duration;
            sdata.(uci).(part_now).frate = pfrate;
            sdata.(part_now).duration = part_duration;         
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Waveform data        
            if numel(pspx)<1 % if there are no spikes there is no point in continuing
                continue
            end

            % get the waveforms for this cluster         
            sdata.(uci).(part_now).waveform_max = [NaN NaN NaN NaN];
            sdata.(uci).(part_now).waveform_width = [NaN NaN NaN NaN];            
            if config.waveform_analyses            
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

                    sdata.(uci).(part_now).waveform_mean{w} = ch;
                    sdata.(uci).(part_now).waveform_stdv{w} = chs;
                    sdata.(uci).(part_now).waveform_max(w) = maxval;
                    sdata.(uci).(part_now).waveform_min(w) = postminval;
                    sdata.(uci).(part_now).waveform_width(w) = width;
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
                sdata.(uci).(part_now).channel_snr = rmss ./ rmsn; % accumulate data
                clear waves2 % clear all the waveform data, we don't need it again and it is quite large
            end
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate ratemap and dwellmap   
            sdata.(uci).(part_now).grid_score = NaN;
            sdata.(uci).(part_now).spatial_measures.spatial_information = NaN;
            if config.field_analyses || config.grid_analyses || config.overdispersion_analyses
                % if a previously computed dwellmap exists, use this to save time
                if isfield(sdata.(part_now),'dwellmap')
                    config.dwellmap = sdata.(part_now).dwellmap;
                else
                    config.dwellmap = [];
                end
                [ratemap,dwellmap,~,config,mapdata] = mapDATA([ppox ppoy],[pspx pspy],config);

                % accumulate data
                sdata.(part_now).dwellmap = dwellmap; % add data to structure            
                sdata.(uci).(part_now).ratemap = ratemap; % add data to structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ratemap analyses      
                % ratemap analysis
                spatial_measures = spatialMETRICS(ratemap,dwellmap);
                sdata.(uci).(part_now).spatial_measures = spatial_measures; % add data to structure          

                % place field analysis
                [fieldd,nfields] = getPFIELDS(ratemap,config.frcut,config.minfr,config.arcut);

                % accumulate data
                sdata.(uci).(part_now).field_count = nfields; % add data to structure
                sdata.(uci).(part_now).field_data = fieldd; % add data to structure   
            
                % analyse grid characteristics
                if config.grid_analyses            
                    % create autocorrelation
                    automap = ndautoCORR(ratemap,ratemap,50);

                    % autocorrelation analysis
                    [grid_score,gdata] = gridSCORE(automap);

                    % accumulate data
                    sdata.(uci).(part_now).grid_autocorrelation = automap; % add data to structure
                    sdata.(uci).(part_now).grid_score = grid_score; % add data to structure
                    sdata.(uci).(part_now).grid_metrics = gdata;                      
                end            
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overdispersion analysis     
            if config.overdispersion_analyses
% Taken from: Fenton et al. (2010) Attention-like modulation of hippocampus place cell discharge
                OD_bin_size = 5; % time window of overdispersion analysis

% In the first method, the entire session was divided into 5-sec intervals. For each interval we calculated the expected number of spikes, exp, as:
% [equation in paper]
% where ri is the time-averaged rate at location i, and ti is the time spent in location i during the pass. Only intervals during which exp ? 5.0 AP were 
% used to calculate overdispersion since the overall firing rate of place cells is ~1.0 AP/sec.
                % convert position coordinates to map coordinates and extract firing rate values from map
                poxmap = mapdata.poxnew;
                poymap = mapdata.poynew;
                exp_frate = ratemap(sub2ind(size(ratemap),round(poymap),round(poxmap)));

                % bin spikes into OD_bin_size second windows
                edg = min(ppot):OD_bin_size:max(ppot);
                [~,~,bindex] = histcounts(ppot,edg);

                % bin observed spikes into the same windows
                [obs_spikes,~,~] = histcounts(pspt,edg);
                bindex(bindex==0) = nanmax(bindex(:))+1; % some spikes fall outside our maximum bin, we just ignore these

                % calculate the expected number of spikes for each time window, based on the firing rate map
                exp_spikes = accumarray(bindex,exp_frate.*(1/50))';     
                exp_spikes(end) = []; % remove the maximum bin we made before

% For each selected 5-sec interval, we then calculated z, the normalized standard deviation of obs, the observed number of spikes as:
% [equation in paper]
                % calculate z, the overdispersion metric
                z = (obs_spikes - exp_spikes) ./ sqrt(exp_spikes);
                z(exp_spikes<5) = [];
                overd = nanstd(z)^2;
                overr = corr(obs_spikes(:),exp_spikes(:),'type','Spearman','rows','pairwise');
% z measures the deviation of observed discharge from expected in standard deviation units. Overdispersion in turn is the variance of 
% the z distribution for a set of passes. The outcome of this calculation was found to be indistinguishable from the somewhat different 
% method previously used (Fenton and Muller, 1998).
% z should have values somewhere between 2 and 5, there is no real upper bound though

                % accumulate data
                sdata.(uci).(part_now).over_dispersion_z = z; % add data to structure            
                sdata.(uci).(part_now).over_dispersion = overd; % add data to structure            
                sdata.(uci).(part_now).over_dispersion_r = overr; % add data to structure 
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HD analyses      
            sdata.(uci).(part_now).hd_rayleigh = NaN;
            sdata.(uci).(part_now).hd_density_rayleigh = NaN;
            if config.hd_analyses
                % Head direction analysis
                ai = linspace(0,2*pi,config.hd_bins)'; % angles for binning
                ai = ai(:);
                hd_c = deg2rad(psph); % cell head direction in radians

                if strcmp(config.hd_type,'density')
                    % HD density estimate
                    if isfield(sdata.(part_now),'hd_dwellmap')
                        hd1 = sdata.(part_now).hd_dwellmap;
                    else
                        hd_s = deg2rad(ppoh); % session head direction in radians
                        [hd1] = circ_ksdensity(hd_s,ai,[],config.hd_sigma); % the session head direction       
                        hd1 = hd1 .* pos_tb; % convert session HD to time 
                        sdata.(part_now).hd_dwellmap = hd1;
                    end

                    [hd2] = circ_ksdensity(hd_c,ai,[],config.hd_sigma); % the cell's head direction
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
                    sdata.(uci).(part_now).hd_density_session = hd1; % add data to structure
                    sdata.(uci).(part_now).hd_density_cell = hd3; % add data to structure
                    sdata.(uci).(part_now).hd_density_frate = hd_c; % add data to structure
                    sdata.(uci).(part_now).hd_density_rayleigh = rayleigh; % add data to structure
                    sdata.(uci).(part_now).hd_density_maximum = mx2; % add data to structure                
                    sdata.(uci).(part_now).hd_density_mean = mn2; % add data to structure                
                    sdata.(uci).(part_now).hd_density_stdev = sd2; % add data to structure     

                elseif strcmp(config.hd_type,'histogram')
                    % HD histogram
                    if isfield(sdata.(part_now),'hd_dwellmap')
                        hd1 = sdata.(part_now).hd_dwellmap;
                    else
                        hd1 = hist(deg2rad(ppoh),config.hd_bins); % the session head direction   
                        fh = fspecial('average',[1 config.hd_boxcar]); % boxcar filter                
                        hd1 = imfilter(hd1,fh,'circular','same');            
                        hd1 = hd1 .* pos_tb; % convert session HD to time
                        sdata.(part_now).hd_dwellmap = hd1;
                    end                

                    hd2 = hist(deg2rad(psph),config.hd_bins); % the cell's head direction
                    hd2 = imfilter(hd2,fh,'circular','same');            
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
                    sdata.(uci).(part_now).hd_session = hd1; % add data to structure
                    sdata.(uci).(part_now).hd_cell = hd3; % add data to structure
                    sdata.(uci).(part_now).hd_frate = hd_c; % add data to structure
                    sdata.(uci).(part_now).hd_rayleigh = rayleigh; % add data to structure
                    sdata.(uci).(part_now).hd_maximum = mx2; % add data to structure
                    sdata.(uci).(part_now).hd_mean = mn2; % add data to structure                
                    sdata.(uci).(part_now).hd_stdev = sd2; % add data to structure
                end
            end
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike phase analyses   
            if config.spike_phase_analyses
                % Spike phase analysis
                % phase preference plot
                lfpf = lfp_theta;
                sdata.theta = int16(lfpf); % add data to structure    
                t = lfp_t;

                % estimate theta phase, based on lfp data
                [phase_out,~,~] = Phase([t lfpf],double(pspt));
                spp = phase_out(:,2);

                % bin the theta phase data
                ai = deg2rad(-180:5:540); % doing this means we have bins symmetrical around zero
                ai = ai(:);
                spp2 = [spp; (spp+2*pi)];
                spp2 = spp2(:);
                yi = histc(spp2,ai);
                phprobs = circ_ksdensity(spp2,ai);

                % directional analyses on phase data
                rayleigh2r = circ_r(spp(:)); % rayleigh vector length
                rayleigh2p = circ_rtest(spp(:)); % rayleigh test for non-uniformity
                mx2p = ai(yi == max(yi)); % preferred angle (location of max frate)

                % example sine wave
                swav = cos(ai);
                swav = ((swav ./ max(swav))+abs(min((swav ./ max(swav))))).*max(yi);
                phmu = circ_mean(spp);

                % accumulate data
                sdata.(uci).(part_now).phase_mean = phmu; % add data to structure
                sdata.(uci).(part_now).phase_r = rayleigh2r; % add data to structure
                sdata.(uci).(part_now).phase_p = rayleigh2p; % add data to structure            
                sdata.(uci).(part_now).phase_max = mx2p(1); % add data to structure
                sdata.(uci).(part_now).spike_phase = spp(:); % add data to structure
                sdata.(uci).(part_now).spike_phase_ideal = swav(:); % add data to structure
                sdata.(uci).(part_now).spike_phase_binned = yi(:); % add data to structure            
                sdata.(uci).(part_now).spike_phase_ksdensity = phprobs(:); % add data to structure 
            end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inter-spike interval analyses   
            if config.spike_isi_analyses || config.spike_autocorrelation_analyses
% Burst index also from Mizuseki et al. (2012) Activity Dynamics and Behavioral Correlates of CA3 and CA1 Hippocampal Pyramidal Neurons 
                % inter-spike interval
                [~,idata,isis] = getISIhalfwidth(pspt);

                % burst index
                bindx = zeros(size(pspt));
                bindx([isis; NaN] < 6 | [NaN; isis] < 6) = 1; % bindx is an index of all spikes sharing an isi less than 6ms
                bindx = logical(bindx);
                burst_index = sum(bindx) / numel(pspt);

                sts = regionprops(bindx,'Area');
                burst_length_median = nanmedian([sts.Area].');
                burst_length_mean = nanmean([sts.Area].');

                % accumulate data
                sdata.(uci).(part_now).isi_data = idata;         
                sdata.(uci).(part_now).burst_index = burst_index;
                sdata.(uci).(part_now).burst_length_median = burst_length_median;
                sdata.(uci).(part_now).burst_length_mean = burst_length_mean;
            end
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Theta and bursting analyses   
            if config.spike_autocorrelation_analyses
% Theta index from Mizuseki et al. (2012) Activity Dynamics and Behavioral Correlates of CA3 and CA1 Hippocampal Pyramidal Neurons 
% Burst index also from Mizuseki et al. (2012) Activity Dynamics and Behavioral Correlates of CA3 and CA1 Hippocampal Pyramidal Neurons 
% Theta skipping index from Brandon et al. (2013) Segregation of cortical head direction cell assemblies on alternating theta cycles - theta_skip is bound between -1 and 1, and higher values indicate more theta cycle skipping.                      
                % Spike autocorrelation - theta analaysis           
                % 250ms autocorrelation
                [tms1,corrdata1] = spikeINTERVALS_v2(pspt,250,1);  
                sdata.(uci).(part_now).theta_train = [corrdata1(:) tms1(:)]; % accumulate data

                % 400ms autocorrelation    
                [tms3,corrdata3] = spikeINTERVALS_v2(pspt,400,1);  
                corrdata3(tms3 == 0) = 0;         
                uni_corr = corrdata3 ./ nansum(corrdata3(:)); % we want the sum of probability to be unity

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
                ft = fittype('(a1*(cos(w*x)+1) + a2*(cos(0.5*w*x)+1) + b) * (exp(-abs(x)/t1)) + (c.*exp(-x^2/(t2^2)))','independent','x','dependent','y');
                opts = fitoptions('Method','NonlinearLeastSquares');
                opts.Algorithm = 'Trust-Region';
                opts.Display = 'Off';
                m = nanmax(prob_skip);
                opts.Lower = [0 0 0 -m 0 0 pi*10];
                opts.Upper = [m m m m 5 0.05 pi*18];
                opts.StartPoint = [0.865 0.042 0.931 0.674 0.567 0.226 0.988].*m;
                [fitresult,~] = fit(xData,yData,ft,opts);
                p1x = (2*pi)/fitresult.w;
                p1 = fitresult(p1x);
                p2x = (4*pi)/fitresult.w;
                p2 = fitresult(p2x);
                theta_skip = (p2-p1) ./ max([p2 p1]);                 

                % accumulate data        
                sdata.(uci).(part_now).theta_train_long = [corrdata3(:) tms3(:) uni_corr_filt(:)];
                sdata.(uci).(part_now).theta_index = theta_index;
                sdata.(uci).(part_now).theta_ratio = theta_ratio;
                sdata.(uci).(part_now).theta_skip = theta_skip;
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speed modulation analyses   
            if config.speed_analyses
%     Taken from Kropff, Carmichael, Moser and Moser (2015) Speed cells in the medial entorhinal cortex
%     Instantaneous firing rate was obtained by dividing the whole session into 20-ms bins, coinciding 
%     with the frames of the tracking camera. A temporal histogram of spiking was then obtained, smoothed 
%     with a 250-ms-wide Gaussian filter.
                [ppotu,~] = unique(ppot); % sometimes dacqUSB repeats time stamps at the end of a file
                edg = ppotu(1)-.01 : .02 : ppotu(end)+.01; % 20ms time bins spanning the session
                shist = histcounts(pspt,edg); % spike histogram
                fh = fspecial('gaussian',[1 13]); % gaussian filter, they say 250ms wide, but samples are every 20ms, so I have used 13 samples or 260ms instead
                shist = imfilter(shist,fh,'replicate','same'); % smooth histogram
                shist = shist ./ 0.02; % convert spikes to firing rate

                % here we bin the speed data and calculate the mean firing rate per 2 cm/s bin
                ppov = fillmissing(ppov,'nearest'); % at least one speed sample can never be computed, so we need to replace that NaN value
                povidx = fix(ppov ./ 2)+1; % bin the speed data in 2cm/s increments
                povidx = interp1(ppot,povidx,histcents(edg),'nearest');
                A = accumarray(povidx(:),shist(:),[],@nanmean);
                speed_time = accumarray(povidx(:),ones(size(povidx))) .* (1/50); % the total time spent moving at each speed

%     The speed score for each cell was defined as the Pearson product-moment correlation between the cell’s 
%     instantaneous firing rate and the rat’s instantaneous running speed, on a scale from ?1 to 1.
                svals = ((1:length(A)).*2)' - 1; % values of speed bins in cm/s
                fvals = A(:); % mean firing rate for each bin
                sscore = corr(svals,fvals,'type','Pearson','rows','pairwise'); % speed score is the correlation between instantaneous firing rate and speed

%     Because of the variability in baseline and slope, a simple or normalized average of speed cell activity would not properly 
%     capture the population behaviour. To obtain a better normalization method, we applied to any firing rate measure f of a 
%     speed cell expressed in Hz the linear transformation [equation in paper]
%     where A (Hz) and B (cm?1) are the y intercept and slope of the cell’s speed 
%     tuning. The 50 value is given in cm s?1. This linear transformation aims to achieve for every cell a normalized dimensionless 
%     activity of 0 when the rat is still and 1 when it runs at 50 cm s?1, allowing for proper population averaging.
                P = polyfit(svals,fvals,1);
                fvals2 = (fvals - P(2)) / (P(1) * 50);

                % accumulate data
                sdata.(uci).(part_now).speed_time = speed_time;
                sdata.(uci).(part_now).speed_score = sscore;
                sdata.(uci).(part_now).speed_slope = P(1);
                sdata.(uci).(part_now).speed_intercept = P(2);
                sdata.(uci).(part_now).speed_curve = [svals fvals];
                sdata.(uci).(part_now).speed_curve_normalised = [svals fvals2];
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Refractory period analyses      
            if config.spike_autocorrelation_analyses
                % Spike autocorrelation - refractory period analysis
                [tms2,corrdata2] = spikeINTERVALS_v2(pspt,25,0.5);     

                % refractory period violation analyses
                tmt = tau_r - tau_c; % max time in which we can observe a noise spike
                nspikes = numel(ppox); % number of cluster spikes
                frate = nspikes ./ part_duration; % firing rate of cluster in Hz
                nrpv = sum(isis <= tau_r); % number of refractory period violations
                prpv = nrpv ./ nspikes; % refractory period violations as a proportion of all spikes

                % Method 1 taken from Dan Hill's UltraMegaSort
                % Spikes occurring under the minimum ISI are by definition false positives, their total rate relative to the rate of the neuron overall gives the estimated false positive rate.
                vtime = 2 * nspikes * tmt; % total time available for violations - there is an available window of length (tau_r - tau_c) after each spike.
                vrate = nrpv/vtime;
                fp_rate1 = vrate/frate;
                if fp_rate1 > 1
                    fp_rate1 = 1; % happens sometimes and is actually an error in the assumptions of the analysis
                end

                % Method 2 using equation from Hill et al. (2011) Quality metrics to accompany spike sorting of extracellular signals
                % This method might return imaginary results for large numbers of violations
                k = 2 * tmt * (nspikes^2);
                rts = roots([k -k nrpv*part_duration]);
                fp_rate2 = min(rts);

                % accumulate data
                sdata.tau_r = tau_r;
                sdata.tau_c = tau_c;            
                sdata.(uci).(part_now).refractory_period = [corrdata2(:) tms2(:)];            
                sdata.(uci).(part_now).rpv_total = nrpv;       
                sdata.(uci).(part_now).rpv_proportion = prpv;                   
                sdata.(uci).(part_now).rpv_false_positive1 = fp_rate1;
                sdata.(uci).(part_now).rpv_false_positive2 = fp_rate2;
                sdata.(uci).(part_now).rpv_censored = (tau_c/1000) * spik_count / part_duration; % estimate the fraction of spikes not detected because of the system lockout             
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell type/identity
            if config.cell_analyses
                [ctype,ctypen,ctypeb] = getCELLTYPE(sdata,uci,part_now);
                sdata.(uci).(part_now).cell_type = ctype; % add data to structure                
                sdata.(uci).(part_now).cell_type_num = ctypen; % add data to structure   
                sdata.(uci).(part_now).cell_type_bin = ctypeb; % add data to structure    
            end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures
            %% Overall summary figure
%             assignin('base','sdata',sdata);
%             assignin('base','uci',uci);
%             assignin('base','part_now',part_now);
%             assignin('base','fig_vis','on');
%             assignin('base','save_fig',save_fig);
%             return

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
        
        dispPROGRESS(cc,1:length(clus)); % display progress percentage     
    end

    %% Cluster cross-correlation figure
    if run_figCROSS
        if sum(clus ~= 0) > 1
            [~,~,~] = mkdir([pwd '\klustest\' sdata.combined_name '\figures\']);
            figfile = [pwd '\klustest\' sdata.combined_name '\figures\' sdata.combined_name '_E' num2str(tet) '_cross_correlogram.png'];            
            figCROSS(sdata,tet,figfile,fig_vis); % overall figure with ratemaps etc
        end
    end

    %% Cluster space figure
    if run_figCSPACE
        [~,~,~] = mkdir([pwd '\klustest\' sdata.combined_name '\figures\']);
        figfile = [pwd '\klustest\' sdata.combined_name '\figures\' sdata.combined_name '_E' num2str(tet) '_cluster_space.png'];
        figCSPACE(sdata.(tetstr).fetdata,figfile,fig_vis); % overall figure with ratemaps etc
    end
end

save([pwd '\klustest\' cname '\' cname '_sdata.mat'],'-struct','sdata'); % save session data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save figures if requested
% disp(sprintf('Saving figures...'))
% 
% % get cell names
% start_string = ['r' sdata.rat_num '_' sdata.date]; % the start of every cluster fieldname
% fdnames = fieldnames(sdata);
% cindx = strfind(fdnames,start_string);
% cindx(cellfun('isempty',cindx)) = {NaN}; % fill empty index values with NaN
% cindx = find(cell2mat(cindx)==1); % find where the non-NaNs are
% cell_names = fdnames(cindx);
% 
% % save figures in parallel
% parfor kk = 1:length(cell_names)
%     uci = cell_names{kk};   
%     sst = strsplit(uci,'_'); 
%     clu = str2double(sst{4}(regexp(sst{4},'\d')));   
%     if ~clu
%         continue
%     end
%     
%     for pp = 1:nparts % for every partition
%         if run_figCLUST
%             figCLUST(sdata,cell_names{kk},part_names{pp},fig_vis,save_fig); % overall figure with ratemaps etc
%         end
%     end
% end

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





























