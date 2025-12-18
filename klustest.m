function klustest(varargin)
% klustest analyse cluster properties
%   _  _ _    _  _ ____ ___ ____ ____ ___ 
%   |_/  |    |  | [__   |  |___ [__   |  
%   | \_ |___ |__| ___]  |  |___ ___]  |  
%
%  klustest analyse clustered data and export analysis figure(s)
%
%  klustest(varargin) uses additional settings specified in name 
%  value pairs (see below)
%
%  name-value input options include:
%
%  'tetrodes'      -   Numeric vector that specified the tetrodes to run on 
%                       (i.e. [1 2 3 4 5 6] would run on tetrodes 1 to 6)
% 
%                       Default value is 1:16
%
%  'clusters'      -   Numeric vector that specifies the clusters to run on
%                       Set to 0 to run on all clusters 
%
%                       Default value is 0
%
%  'rname'         -   String or charcter vector specifying the rat name
%                       This is used in the sdata structure so its important 
%                       to ensure this is correct
%
%                       Default value is the name of the parent directory
% 
%  'outname'       -   String or charcter vector specifying the output
%                       filename used by klustest - i.e. figures will be saved
%                       in a folder named this.
%
%                       Default value is 'klustest'.
% 
%  'cname'         -   String or charcter vector specifying the output
%                       filename used by kwiktint - or the name used for
%                       klustakwiked files
%
%                       Default value is 'kwiktint'.
%
%  Notes
%  -----
%  1. If cells are recorded in the same session their spikemaps are 
%      expected to differ while the position data and dwellmap will
%      remain the same. To save time this function can accept a precomputed 
%      output, skipping that computation (which is usually more time 
%      consuming than the spikemap). Simply provide the 'speedlift' output 
%      from rate_mapper to the next iteration of rate_mapper. For some 
%      methods a dwellmap is not computed or would not speed up 
%      computation if provided, so instead a matrix or similar form of 
%      data output is used instead. Thus in some cases (i.e. 'histogram' 
%      method) the speedlift matrix will be identical to the dwellmap, 
%      but in other cases (i.e. 'kadaptive') it will take a different 
%      form and contain different data.
%
%  2. For the 'histogram' method and smoothing method 1, smoothing is 
%      achieved using imgaussfilt. The FilterDomain is set to 'spatial'
%      to ensure convolution in the spatial domain. 'Padding' is set to 
%      a scalar value of 0 as there should be no position or spike data
%      in bins outside the map limits.
% 
%  3. For the 'histogram' method and smoothing method 1, smoothing is 
%      achieved using imgaussfilt. The FilterDomain is set to 'spatial'
%      to ensure convolution in the spatial domain. 'Padding' is set to 
%      a scalar value of 0 as there should be no position or spike data
%      in bins outside the map limits.
%
%
%  Example
%  ---------
%
%  % in a directory with kwiktint outputs, run function using default values
%  klustest()
% 
%  % run function using default values, but only on tetrodes 1 and 5
%  klustest('tetrodes',[1 5])
% 
%  % run function using default values, all specified
%  klustest('tetrodes',1:16,'clusters',0,'rname','RAT1','cname','kwiktint','outname','klustest')
%
%  See also kwiktint

% HISTORY:
% version 01.0.0, Release 23/02/16 YTcluanalysisROD majorly revised
% version 01.0.1, Release 24/02/16 added functionality to find session details from .set file
% version 01.0.2, Release 30/03/16 fixed some problems with .set file data because I updated it to handle multiple sessions
% version 01.0.3, Release 05/08/16 renamed, fixed some errors
% version 01.1.0, Release 08/08/16 modified readDACQDATA to output waveforms and added post processing of this info
% version 01.2.0, Release 08/08/16 readDACQDATA working, added postprocessing of mtint
% version 01.3.0, Release 09/08/16 main plots working
% version 01.4.0, Release 10/08/16 options added for different styles of plot
% version 02.0.0, Release 11/08/16 added theta autocorrelation sine estimate (O'Mara 2015)
% version 02.1.0, Release 12/08/16 added refractory period violations (Navratilova and McNaughton (2016)
% version 02.1.1, Release 13/08/16 added saving data, figures
% version 02.2.0, Release 14/08/16 added cell type identification
% version 03.0.0, Release 15/08/16 cluster quality assessment added in clusterQUALITY
% version 03.1.0, Release 16/08/16 added cluster space figures
% version 03.1.1, Release 17/08/16 fixed phase plot and phase map
% version 03.2.0, Release 18/08/16 modified clusterQUALITY to deal better with missing channels
% version 03.3.0, Release 18/08/16 added the option to ignore position data, fixed bug with waveform plot
% version 03.3.1, Release 19/08/16 make it so none of the position figures are made if there is no position data, this should be faster
% version 03.4.0, Release 22/08/16 fixed issues with cluster space plot, concentrate on first feature, added legend
% version 03.5.0, Release 23/08/16 added cluster space subplot to cell figure, for Ele, this uses the first clustering feature, the 1st and 2nd highest amplitude channel
% version 04.0.0, Release 25/06/16 added tetrode input to readalldacqdata, this prevents it trying to open unnecessary files, fixed bug in the detection of dead/lfp channels
% version 04.1.0, Release 25/08/16 fixed a minor bug in histograms when there is only one spike
% version 05.0.0, Release 22/11/16 adapted from KlusterAnalysisTINT
% version 05.1.0, Release 30/03/17 started editing the newer version with more emphasis on sdata structure
% version 06.0.0, Release 30/03/17 adapted to new partitioning system using a part_config cell array, function can now split based on any inputs
% version 07.0.0, Release 31/03/17 started moving figure work into subfunctions to contain them, each uses the sdata structure
% version 07.1.0, Release 01/04/17 finished most of the figure subfunctions
% version 08.0.0, Release 02/04/17 added kluspartfig, which generates a summary figure containing all the parts, and correlations between them
% version 09.0.0, Release 05/04/17 fixed multiple problems in clusterQUALITY, including wrong lratios and reliance on kwiktint files, now using clusterQUALITY_v2
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
% version 16.0.0, Release 04/04/19 moved part_config guts to later in the function, replaced part_config structure to be a table, now saved as .mat
% version 16.0.1, Release 04/04/19 replaced getTRODES and getCNAMES with a combined version called getTRODES
% version 16.0.2, Release 05/04/19 moved pixel ratio handling to postprocessDACQDATA
% version 16.0.3, Release 07/04/19 part indexing is now built logically instead of in a loop (should be faster)
% version 16.0.4, Release 08/04/19 updated spatial mapping analyses
% version 16.1.0, Release 09/04/19 moved HD bining to mapHD, theta curve fitting to getTHETAintrinsic, refractory analysis to getRPVcount, speed analysis to getSPEEDmod      
% version 16.1.1, Release 09/04/19 changed getCELLTYPE for getCELLTYPE2 which outputs a binary approach to cell typing
% version 16.2.0, Release 10/04/19 moved cluster quality to getDACQDATA
% version 16.2.1, Release 10/04/19 better handling of waveforms
% version 16.2.1, Release 13/04/19 figure functions finished, figPART and figCLUST
% version 16.3.0, Release 15/04/19 Finished converting sdata to table format, variables have to be preallocated by subfunction addToTable
% version 16.3.1, Release 15/04/19 progress is now displayed using fileExchange function 'ProgressBar'
% version 16.3.2, Release 15/04/19 improved handling of keypresses in figKEY, added manual definition of start/end keys, part_config now contains files required for each part
% version 16.3.3, Release 16/05/19 progress shown using custom code 'looper'
% version 16.4.0, Release 20/06/19 fixed error when part_config is first created because Matlab's use of cell arrays in tables changed with 2019b
% version 17.0.0, Release 15/04/19 added exception for weird characters in saveKEY, fixed error in figKEYS where timestamps were not concatenated correctly
% version 17.1.0, Release 15/04/19 overhaul of digital input processing, added manageKEYS which can be used via prepKEYS without klustest
% version 17.2.0, Release 17/04/19 added name value pair arguments
% version 18.0.0, Release 13/02/23 major overhaul started to accept Neuralynx formats, merging different versions of klustest, improving figures and speed, improving sdata
% version 18.1.0, Release 13/02/23 data loading now done by get_tets_for_klustest, get_pos_for_klustest, get_clu_for_klustest, etc, this allows multiple data formats to be loaded
% version 18.2.0, Release 13/02/23 added ability to load 3D reconstructed position data, simplified 3D HD calculation
% version 18.3.0, Release 13/02/23 improved part_config and part figure which is now made by overhauled klustfig_part
% version 18.3.1, Release 14/02/23 improved comments and command line messaging
% version 18.4.0, Release 14/02/23 replaced mapping function with rate_mapper
% version 18.4.1, Release 14/02/23 fixed bug in theta phase analysis
% version 18.5.0, Release 14/02/23 removed bloat from sdata, removed unnecessary preallocation, pdata is now stored in sdata as a custom property structure
% version 18.5.1, Release 15/02/23 fixed bugs in klustfig_part, function seems to run smoothly in cases where there are no spikes now
% version 18.6.0, Release 15/02/23 added analysis_log so analysis steps can be recorded and easily checked
% version 18.6.1, Release 18/02/23 added limits to klustfig_part so it doesn't plot every waveform or point in the feature space
% version 18.6.2, Release 19/02/23 added position data to first cluster of a session, added 3D HD to get_pos_for_klustest('reconstruction')
% version 19.0.0, Release 22/06/25 added support for Axona data formats 
% version 19.0.1, Release 23/06/25 improved figures and comments
% version 19.0.2, Release 24/06/25 simplified sdata and added bdata tables
% version 19.0.3, Release 24/06/25 simplified autocorrelations and sdata contents
% version 20.0.0, Release 29/11/25 Updates for GitHub release
% version 20.0.1, Release 29/11/25 datesrt depreciated, replaced with datetime
% version 20.0.2, Release 16/12/25 updated get_iso_for_klustest to use logarithmic scale for histograms
% version 20.0.3, Release 16/12/25 removed inset CDFs in isolation graph in klustfig_part
% version 20.0.3, Release 16/12/25 updated filenames for cross platform flexibility
% version 20.0.4, Release 17/12/25 fixed bug where incomplete cluster quality could be loaded
% version 21.0.0, Release 17/12/25 readtable/writetable no longer viable for part_config, changing this to a structure saved in a .json makes the code no longer backward compatible
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2025 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    p = inputParser;
    addParameter(p,'tetrodes',1:16,@(x) ~isempty(x) && ~all(isnan(x(:))) && isnumeric(x));  
    addParameter(p,'clusters',0,@(x) ~isempty(x) && ~all(isnan(x(:))) && isnumeric(x)); 
    addParameter(p,'screening',0,@(x) ~isscalar(x) && isnumeric(x));     
    addParameter(p,'cname','kwiktint',@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    addParameter(p,'outname','klustest',@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x));     
    dlist = strsplit(pwd,filesep);
    rname = dlist{end-1};
    addParameter(p,'rname',rname,@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    parse(p,varargin{:});   
    config = p.Results;  

    %% overrides
    override.part_config = 0; % 0 = use precomputed part_config file when available
    % override.use_groups = 1;
    override.quickstart = 0; % 0 = use precomputed LFP and cluster quality when available
    config.spatial_shuffles = 1; % 1 = calculate probability of grid score, spatial info, rayleigh vector
    config.skipfigs = 0; % 1 = skip making a figure if it exists already    
    config.save_figs = 1; % 1 = make and save figures (VERY time consuming, takes about 90% of klustest's time)
    config.fast_figures = 1; % 1 = 4-5x faster figure saving, but at a slightly lower quality
    config.fig_res = 250; % resolution at which to save figures
    config.n_waveforms = 250; % number of waveforms to show
    config.max_plot_spikes = 1e5; % set to 0 to show all spikes in spike & position plot
    config.show_figs = 'off'; % display figures before saving
    
    %% Map settings
    mapset.ppm          = 300; % pixels per meter, can be scalar (same value usd for all sessions) or Nx1, were N = number of sesions
    mapset.jumpcut      = 5; % standard deviations, jumps in the position data with a zscored distance between them greater than this will be removed and interpolated
    mapset.jumpwindow   = 10; % number of samples over which to smooth jump detection, higher values will remove more data before/after jumps, increase this if outliers are being left behind
    mapset.method       = 'histogram';
    mapset.binsize      = 25; % (mm) firing rate map bin size
    mapset.ssigma       = 30; % (mm) firing rate map smoothing sigma
    mapset.padding      = 10; % (mm) how much to pad the edges of the firing rate map
    mapset.mindwell     = 0.01; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    mapset.mindist      = 40; % (mm, default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    mapset.smethod      = 2; % smoothing method, 1 = before division, 2 = after, 3 = no smoothing 
    mapset.zcut         = 1; % z-score for field detection in z-scored firing rate maps   
    mapset.frcut        = 1; % place fields must have a peak firing rate at least greater than this value
    mapset.arcut        = 400; % cm2, place fields must have a total area greater than this value
    mapset.fix_aspect   = 1; % if set to 1, position plots will be rotated so that they are always landscape (the orientation of the figures and thus the best for visualisation)
    mapset.wave_window  = [-0.2 0.8]; % (default [-0.25 0.75], Axona = [-0.2 0.8]) ms, the time window over which to load and plot waveforms (only used when loading Phy data, but used for all plotting)
    
    % HD settings - settings to use when generating head direction maps
    mapset.hd_displace  = 0; % (default 0) if set to 1, will use displacement (movement direction) instead of head direction
    mapset.hd_type      = 'histogram'; % (default 'histogram') enter 'density' for a circular kolmogorov smirnov density estimate plot, enter 'histogram' for the traditional histogram polar plot
    mapset.hd_bins      = 60; % (default 64) the number of bins to use when computing HD plot
    mapset.hd_sigma     = 0.04; % (default 2) the standard deviation of the gaussian kernel used in HD circular density estimate
    mapset.hd_boxcar    = 3; % (default 3) the number of bins over which to compute the HD histogram boxcar average    
    
    %% Data formats
    % % if cluster cut with kwikcut
    % formats.pos         = 'Neuralynx';
    % formats.clu         = 'Neuralynx';
    % formats.tet         = 'Tint';
    % formats.set         = 'Neuralynx';    
    % formats.spk         = 'Neuralynx';
    % formats.lfp         = 'Neuralynx';
    % formats.iso         = 'klustakwik';
    
    % if cluster cut with kwiktint
    formats.pos         = 'kwiktint';
    formats.clu         = 'kwiktint';
    formats.tet         = 'kwiktint';
    formats.set         = 'kwiktint';    
    formats.spk         = 'kwiktint';
    formats.lfp         = 'kwiktint';
    formats.iso         = 'kwiktint';   
    formats.pos         = 'kwiktint';
    formats.front_led_color = 1; % 1 = red, 2 = green, 3 = blue
    formats.back_led_color = 2; % 1 = red, 2 = green, 3 = blue
    formats.led_angle_offset = 0; % CCW offset of LEDs on the head

%%%%%%%%%%%%%%%% PARTITION SETTINGS
    % This function is assumed to follow spike sorting using kwiktint. Kwiktint
    % can merge multiple recording sessions into one output. Klustest
    % recognises this can and redivide the data back up into indiviual sessions
    % or 'parts' for analysis and visualisation.
    % The next section deals with how this division or partitioning works.

    % Each part is defined by a set of characteristics, they each need a name,
    % which is used when saving the output. Because these are used as structure
    % field names they cannot start with a number.
    % Common names might be 'session1' or 'part_3' or 'Arena_1'

    % They need a method for partitioning - do we want to split the merged file
    % up according to recordings? Do we want to keep some recordings combined?:
    % 1 = combine everything
    % 2 = take a recording session
    % 3 = use digital inputs

    % If the method above is 2 or 3 we also need to specify which intervals
    % to use as the start and end of the trials we want to extract, or, more
    % commonly, which recording sessions should be combined into a part:
    % if method = 1 this value is ignored
    % if method = 2 this should specify which recording session(s) to include
    % if method = 3 this should specify which digital input pairs to include (inf = use all) 

    % Lastly, if the method is set to 3 (digital inputs), the intervals list
    % tells us which digital input key pairs to use for the start/end of the trials
    % we want to include, but we still need to specify what those digital input key
    % pairs are, some people use 's' and 'e' to mark the start and end of each
    % trial for example. If we are not using method 3, this value is ignored.
    % keypresses used to delineate trial starts and ends {start,end} can be integers 
    % or characters i.e. {1,2} or {'s','e'}, {'1','2'} is the same as {1,2}

    % for example, here we create a part named 'session1', it will take data
    % from a recording session (method 2) which is the first recording made
    % (interval 1) and we don't need any interval keys.
    % part_name                               = 'session1';
    % part_config.(part_name).method          = 2;
    % part_config.(part_name).intervals       = [1];
    % part_config.(part_name).interval_keys   = {};

    % Here we create a part named 'session2', it will take data
    % from recording sessions (method 2) which are the second and third
    % recordings made (intervals 2 & 3) and we don't need any interval keys.
    % part_name                               = 'session1';
    % part_config.(part_name).method          = 2;
    % part_config.(part_name).intervals       = [2 3];
    % part_config.(part_name).interval_keys   = {};

    % Here we create a part named 'session3', it will take data
    % from any recording session that fall between keypresses (method 3)
    % The interval keys it will use to specify the start and end of each trial
    % are 's' for starts and 'e' for ends, it will combine the trials
    % 1,2,3,4,5,6,10,11 and 15:
    % part_name                               = 'session3';
    % part_config.(part_name).method          = 3;
    % part_config.(part_name).intervals       = [1 2 3 4 5 6 10 11 15];
    % part_config.(part_name).interval_keys   = {'s','e'};

    % space for actual part specification:
    part_config = struct;
    part_name                               = 'session1';
    part_config.(part_name).method          = 2;
    part_config.(part_name).intervals       = [1];
    part_config.(part_name).interval_keys   = {};
    % part_name                               = 'session2';
    % part_config.(part_name).method          = 2;
    % part_config.(part_name).intervals       = [2 3 4 5 6];
    % part_config.(part_name).interval_keys   = {};
    % part_name                               = 'session3';
    % part_config.(part_name).method          = 2;
    % part_config.(part_name).intervals       = [7];
    % part_config.(part_name).interval_keys   = {};

    % add additonal information
    pnames = fieldnames(part_config);
    for ff = 1:length(pnames)
        if isscalar(mapset.ppm) 
            part_config.(pnames{ff}).pratio = mapset.ppm;
        else
            part_config.(pnames{ff}).pratio = mapset.ppm(ff);
        end
    end
    
    % save the part_config so that in the future we can just run the function
    % with the same settings 
    jsonText = jsonencode(part_config,'PrettyPrint',true);
    [~,~,~] = mkdir(fullfile(pwd,config.outname)); % create a folder to hold outputs        
    part_config_fname = fullfile(pwd,config.outname,'part_config.json');
    if override.part_config || ~exist(part_config_fname,'file')
        fid = fopen(part_config_fname,'w');
        fwrite(fid,jsonText);
        fclose(fid);
    elseif override.part_config && exist(part_config_fname,'file')
        part_config_fname2 = fullfile(pwd,config.outname,['part_config_' datetime("now",'format','yyyyMMddHHmmss') '.json']);     
        [~,~,~] = movefile(part_config_fname,part_config_fname2);
        fid = fopen(part_config_fname,'w');
        fwrite(fid,jsonText);
        fclose(fid);
    end
    part_config = jsondecode(fileread(part_config_fname));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE DATA
%%%%%%%%%%%%%%%% Tetrodes and sessions
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'); tic;
    t = datetime('now');
    s = string(t,'yyyy-MM-dd HH:mm:ss');    
    disp(sprintf('Running %s at %s...','klustest',s))

    %% find which tetrodes are available
    % When sessions are cluster cut together (which you should do with multiple sessions recording the same cells)
    % kilocut saves the filenames of the individual sessions in the kilo.mat file as well as the tetrodes analysed
    disp(sprintf('Assessing data...'))
    [tets,snames,data_dirs] = get_tets_for_klustest(formats.tet,config);
    config.snames = snames;
    disp(sprintf('\t...%d sessions',size(snames,1)));   
    disp(sprintf('\t...will analyse electrodes: %s',mat2str(tets(:,1)')))  
            
    % get additional session info
    disp(sprintf('\t...session info'));
    [hdata,fname] = get_set_for_klustest(formats.set,config);
    disp(sprintf('\t\t...read: %s',fname));            
    disp(sprintf('\t\t...recording date: %s',hdata.date));        

    % Start to prepare the pdata (part or session data) table
    % To save space and make things easier for the user I have divided the data into two main arrays
    % the pdata structure contains data which applies to the whole session, like positions or dwell maps
    % the sdata table contains data for each cluster/part (one per row). The idea being that tables are easier to concatenate
    % to build a full data set
    pdata = struct;
    pdata.rat = config.rname;
    pdata.date = hdata.date;
    pdata.analysed = s;
    pdata.directory = pwd;
    pdata.tetrodes = tets; % list of tetrodes analysed
    pdata.sessions = size(snames,1); % the number of recording sessions used
    pdata.data_dirs = data_dirs;  
    pdata.snames = snames;
    pdata.mapset = mapset;
    pdata.cname = config.cname;
    pdata.outname = config.outname;

%%%%%%%%%%%%%%%% Position data
    disp(sprintf('\t...positions'));
    pos_srate = 50; % desired sampling rate of position data (Hz)
    [pos,data_intervals,tstart] = get_pos_for_klustest(formats,data_dirs,snames,pos_srate,mapset); % directories
    % pos = [x,y,t,v,hd,ahv] for all data
    % data_intervals = [tstart,tend] one row per recording
    % tstart = the original start time of the first recording
    pdata.pos = pos;
    pdata.pos_srate = pos_srate;
    pdata.tstart = tstart;

%%%%%%%%%%%%%%%% Cluster data
    disp(sprintf('\t...clusters'))
    clus = get_clu_for_klustest(formats.clu,config,tets,data_dirs);
    pdata.clusters = clus;
    
%%%%%%%%%%%%%%%% Spike data
    disp(sprintf('\t...spikes and waveforms'))
    [spk,wav,spk_srate,wavtime] = get_spk_for_klustest(formats.spk,data_dirs,tets,tstart,clus,mapset.wave_window);
    pdata.spike_times = spk;
    pdata.wavtime = wavtime;

%%%%%%%%%%%%%%%% Cluster quality
    disp(sprintf('\t...cluster quality'))
    mname = fullfile(pwd,pdata.outname,'klustest_quickstart.mat');
    run_get_iso = 1;
    if exist(mname,'file') && ~override.quickstart && 1
        matObj = matfile(mname);
        vars = who(matObj);
        if ~any(strcmp(vars,{'isods'})) % if cluster quality is missing from the .mat file
            run_get_iso = 1; % extract it below
        else
            disp(sprintf('\t\t...loading from file'))            
            load(mname,'isods','fets','quals','-mat');

            % if the user previously ran klustest with an incomplete number of
            % tetrodes or a different number of clusters, we need to recompute it
            if numel(isods)==size(tets,1)
                % isods has the correct number of electrodes but does it have the
                % correct number of clusters?
                n_clus = cell2mat(cellfun(@(x) numel(unique(x(x>0))),clus,'UniformOutput',false));
                n_clus_isods = cell2mat(cellfun(@(x) size(x,1),isods,'UniformOutput',false));
                if all(n_clus==n_clus_isods)
                    run_get_iso = 0; % do not extract it below
                else
                    run_get_iso = 1; % extract it below   
                end
            else
                run_get_iso = 1; % extract it below
            end
        end
    end
    if run_get_iso
        [isods,quals,fets] = get_iso_for_klustest(formats.iso,config,tets,clus,data_dirs);
        disp(sprintf('\t\t...saving to file'))                            
        save(mname,'isods','fets','quals','-mat','-v7.3'); 
    end

%%%%%%%%%%%%%%%% LFP data
    disp(sprintf('\t...LFP'))
    mname = fullfile(pwd,pdata.outname,'klustest_quickstart.mat');    
    run_get_lfp = 1;
    if exist(mname,'file') && ~override.quickstart && 1
        matObj = matfile(mname);
        vars = who(matObj);
        if ~any(strcmp(vars,{'lfp'})) % if lfp is missing from the .mat file
            run_get_lfp = 1; % extract it below        
        else
            disp(sprintf('\t\t...loading from file')) 
            load(mname,'lfp','lfp_srate','-mat');    
            run_get_lfp = 0; % do not extract it below  
        end
    end
    if run_get_lfp
        [lfp,lfp_srate] = get_lfp_for_klustest(formats.lfp,data_dirs,tstart,250); % lfp = zscored amplitude, time, theta phase, theta power
        disp(sprintf('\t\t...saving to file'))                            
        save(mname,'lfp','lfp_srate','-mat','-append'); 
    end    

%%%%%%%%%%%%%%%% Part data
    disp(sprintf('\t...part config'))
    part_names = fieldnames(part_config);
    nparts = length(part_names);
    for pp = 1:nparts
        switch part_config.(part_names{pp}).method
            case {1} % if the method is whole session (i.e. everything recorded)
                part_config.(part_names{pp}).interval_times = data_intervals(:);
                % include all the periods where data were recorded

            case {2} % if the method is recording session
                part_config.(part_names{pp}).interval_times = data_intervals(part_config.(part_names{pp}).intervals,:);
                % include all the recording periods listed in part_intervals

            case {3} % if the method is intervals (keypresses)
                % still to add, use manageKEYS and Nlx2MatEV
                keyboard
                
        end
        part_config.(part_names{pp}).part_duration = sum(part_config.(part_names{pp}).interval_times(:,2)-part_config.(part_names{pp}).interval_times(:,1));        
    end    
    pdata.part_config = part_config; % a copy of the part_config, this will be extended to include more data that the saved part_config though, table format
    disp(sprintf('\t\t...done'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERS
    disp(sprintf('Analysing clusters...'))
    sdata = table;
    bdata = table;

    % Kilosort mainly ignores tetrodes and just assigns cluster IDs across 
    % all available channels (even though we well it not to cluster across 
    % tetrodes). So, in this loop we run through every unique cluster returned
    % by kilosort and analyse its activity
    for tt = 1:size(tets,1) % for every tetrode       
        if tets(tt,2)==0 % tetrode
            disp(sprintf('\tTetrode %d | %d of %d (%.f%%)',tets(tt,1),tt,size(tets,1),tt/size(tets,1)*100))
        elseif tets(tt,2)==1 % stereotrode
            disp(sprintf('\tStereotrode %d | %d of %d (%.f%%)',tets(tt,1),tt,size(tets,1),tt/size(tets,1)*100))
        end
        clu = clus{1,tets(tt)}; % clusters on this tetrode
        clusters = unique(clu(clu>0)); % list on nonzero clusters
        spt = pdata.spike_times{1,tets(tt)}; % spike times for this tetrode
        wav_now = wav{1,tets(tt)}; % waveforms for this tetrode
        
        if isempty(clusters)
            disp(sprintf('\t\tno data clusters'))            
        end
        for cc = 1:length(clusters) % for every cluster               
            disp(sprintf('\t\tcluster %d of %d (%.f%%): ',clusters(cc),length(clusters),cc/length(clusters)*100))

            % The next loop focuses on dividing the data into its different parts 
            % and analysing the cluster within each of these. The results are then 
            % added to sdata_temp which will be concatenated with sdata at the end    
            for pp = 1:nparts % for every part               
                sdata_temp = table; % temporary table for this cluster
                bdata_temp = table;
                part_now = part_names{pp};
                if pp==1; disp(sprintf('\b%s ',part_now)); else; disp(sprintf('\b| %s ',part_now)); end

                % cut the position data to include only this part
                interval_times = part_config.(part_now).interval_times;
                if ~isfield(pdata,part_now) || ~isfield(pdata.(part_now),'pox')
                    pos = pdata.pos;
                    pot = pos.pot;
                    pindax = logical(sum(pot' >= interval_times(:,1) & pot' <= interval_times(:,2),1));
                    pdata.(part_now).pot = pos.pot(pindax,1); % pos time for this part            
                    pdata.(part_now).pox = pos.pox(pindax,1); % pos x for this part
                    pdata.(part_now).poy = pos.poy(pindax,1); % pos y for this part
                    pdata.(part_now).poh = pos.poh(pindax,1); % pos (yaw) HD for this part
                    pdata.(part_now).pod = pos.pod(pindax,1); % pos (yaw) HD for this part, estimated using displacement                  
                    pdata.(part_now).pov = pos.pov(pindax,1); % velocity for this part
                    pdata.(part_now).pav = pos.poa(pindax,1); % pos (yaw) AHV for this part
                end
                ppot = pdata.(part_now).pot; % pos time for this part            
                ppox = pdata.(part_now).pox; % pos x for this part
                ppoy = pdata.(part_now).poy; % pos y for this part
                ppoh = pdata.(part_now).poh; % pos (yaw) HD for this part
                ppod = pdata.(part_now).pod; % pos (yaw) HD for this part, estimated using displacement                 
                ppov = pdata.(part_now).pov; % velocity for this part
                ppoa = pdata.(part_now).pav; % pos (yaw) AHV for this part
                    
                % cut the spike data to include only this part and cluster
                sindax = logical(sum(spt' > interval_times(:,1) & spt' < interval_times(:,2) & (clu==clusters(cc))',1)); % spikes                    
                pspt = spt(sindax); % spike time for this part
                sidx = knnsearch(ppot,pspt); % nearest neighbour in position data for every spike
                pspx = ppox(sidx); % spike x for this part
                pspy = ppoy(sidx); % spike y for this part
                psph = ppoh(sidx); % spike HD for this part
                pspd = ppod(sidx); % spike HD for this part, estimated using displacement 
                %pspv = ppov(sidx); % spike velocity for this part (currently not needed)
                pspa = ppoa(sidx); % spike AHV for this part             
                
                % accumulate data in sdata_temp
                sdata_temp.rat = {pdata.rat}; % the rat number/name
                sdata_temp.date = {pdata.date};
                sdata_temp.partn = uint8(pp);
                sdata_temp.dir = [pwd filesep];
                if isempty(bdata) || ~any(bdata.partn==sdata_temp.partn) % this part has not been added to bdata yet
                    bdata_temp = sdata_temp; % will hold behaviour data in table format
                end
                if tets(tt,2)==0 % tetrode
                    uci = ['uci_' rname '_' pdata.date '_t' num2str(tets(tt,1)) '_c' num2str(clusters(cc))]; % unique cell identifier - a string with [rat name, session date, electrode, cluster];
                elseif tets(tt,2)==1 % stereotrode
                    uci = ['uci_' rname '_' pdata.date '_s' num2str(tets(tt,1)) '_c' num2str(clusters(cc))]; % unique cell identifier - a string with [rat name, session date, electrode, cluster];
                end
                sdata_temp.uci = {uci};
                sdata_temp.tetrode = uint8(tets(tt,1)); 
                sdata_temp.electrode = uint8(tets(tt,2));                
                sdata_temp.cluster = uint8(clusters(cc)); % the index of the cluster in Phy format (number assigned to cluster in Phy, may not start at zero, may have gaps etc)
                sdata_temp.spt_pot_index = {uint64(sidx)}; % this index can be used to get the spike values from the position data values i.e. pdata.(part_now).pox(sdata.spike_index)
                sdata_temp.spike_times = {single(pspt)}; % the actual spike times, single should be fine as spike times are sampled at 32kHz
                sdata_temp.nspikes = numel(pspx);
                sdata_temp.frate = numel(pspx) / (numel(ppot)*(1/pos_srate));           
                inow = isods{1,tets(tt)};
                sdata_temp.isod = inow(cc,:); % cluster quality [isolation distance, lratio]
                if isempty(bdata) || ~any(bdata.partn==sdata_temp.partn) % this part has not been added to bdata yet
                    bdata_temp.pos = { single([ppot(:) ppox(:) ppoy(:) ppoh(:) ppov(:) ppoa(:) ppod]) };
                end                               
                
%%%%%%%%%%%%%%%% Waveform data
                nch = size(wav_now,2);
                mxs = NaN(nch,1);
                wav_means = NaN(4,size(wav_now,1));
                wav_stds = NaN(4,size(wav_now,1));  
                wav_now_clus = cell(1,4);
                for ww = 1:nch % for every channel
                    wav_now_clus{ww} = squeeze(wav_now(:,ww,sindax))'; % should be [spikes x samples]
                    wav_means(ww,:) = mean(wav_now_clus{ww},1,'omitnan');
                    wav_stds(ww,:) = std(double(wav_now_clus{ww}),[],1,'omitnan');
                    mxs(ww,1) = max(wav_means(ww,:),[],'omitnan');
                end
                [~,widx] = sort(mxs,'descend','MissingPlacement','last'); % sort from largest > smallest waveform
                max_wav_means = wav_means(widx,:);
                max_wav_stds = wav_stds(widx,:);
                max_wav_mxs = mxs(widx);

                % use the waveform with the largest amplitude to calculate wave characteristics
                wnow = max_wav_means(1,:);
                [~,maxvalt] = max(wnow);
                wnow(1:maxvalt) = NaN;
                [minval,minvalt] = min(wnow,[],'omitnan');
                width_ms = (minvalt - maxvalt) * (1/spk_srate) * 1e03; % waveform width in ms

                % accumulate data in sdata_temp
                sdata_temp.wave_width = single(width_ms);
                sdata_temp.wave_amps = int16(mxs(1));
                sdata_temp.wave_mean = { int16(wav_means) };
                sdata_temp.wave_std = { int16(wav_stds) };
                %pdata.wavtime = (( 1:size(wav_now,1) ) -8 ) * (1/spk_srate) * 1e03; % waveform samples, converted to seconds, then microseconds

%%%%%%%%%%%%%%%% Firing rate map
                % This function incorporates a number of different methods for generating firing rate maps and you can check the function for detailed descriptions of these
                % my favourite is the method proposed by Leutgeb et al. (2005) Independent Codes for Spatial and Episodic Memory in Hippocampal Neuronal Ensembles
                % https://doi.org/10.1126/science.1114037
                % although the fastest method is 'histogram' as this is just a straight binning method followed by smoothing, so it is the default
                speedlift = [];
                if isfield(pdata,part_now)
                    if isfield(pdata.(part_now),'speedlift')
                        speedlift = pdata.(part_now).speedlift;
                    end
                end
                rmset = mapset;
                pos = [ppox ppoy].*10;
                spk = [pspx pspy].*10;
                rmset.maplims = [min(pos(:,1)) min(pos(:,2)) max(pos(:,1)) max(pos(:,2))]; 
                rmset.srate = pos_srate;
                [ratemap,dwellmap,spikemap,~,speedlift] = rate_mapper(pos,spk,rmset,speedlift);         

                % Skaggs spatial information content (bits per second)
                % Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
                % https://onlinelibrary.wiley.com/doi/epdf/10.1002/%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K
                % Markus et al. (1994) Spatial information content and reliability of hippocampal CA1 neurons: effects of visual input
                % https://onlinelibrary.wiley.com/doi/epdf/10.1002/hipo.450040404
                pr = dwellmap ./ sum(dwellmap(:),'omitnan'); % dwell time probability
                ro = sum(ratemap(:) .* pr(:),'omitnan'); % overall firing rate
                si = sum(pr(:) .* (ratemap(:)./ro) .* log2(ratemap(:)./ro),'omitnan'); 

                % Skaggs spatial information content (bits per spike)
                % From Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
                % The information rate given by formula (1) is measured in bits per second. If it is
                % divided by the overall mean ring rate of the cell (expressed in spikes per second),
                % then a different kind of information rate is obtained, in units of bits per spike|let us
                % call it the information per spike. This is a measure of the specificity of the cell: the
                % more grandmotherish" the cell, the more information per spike. 
                sis = si ./ ro;

                % Sparsity
                % From Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
                % The sparsity measure is  an adaptation  to space of  a formula invented by Treves and Rolls (1 99 1); the adaptation measures the 
                % fraction of the environment  in which a cell  is active. Intuitively, a sparsity of, say, 0.1 means that the place field of the cell 
                % occupies 1/10 of the area the rat traverses.
                sp = (sum(pr(:).*ratemap(:),'omitnan').^2) ./ (sum(pr(:).*(ratemap(:).^2),'omitnan'));            

                % Spatial Coherence (Cacucci et al. 2007)
                % https://dx.doi.org/10.1523%2FJNEUROSCI.1704-07.2007
                % The spatial coherence for each firing rate map was computed as the mean correlation between the firing rate of each bin with the 
                % aggregate rate of the 24 nearest bins.
                meanf = ones([5 5]);
                meanf(3,3) = 0;
                U = imfilter(ratemap,meanf,0,'same','conv');
                cohe = corr(ratemap(:),U(:),'rows','pairwise','type','Pearson');            

                % accumulate data
                sdata_temp.ratemap = { single(ratemap) };
                sdata_temp.spikemap = { single(spikemap) };
                sdata_temp.spatial_info = [si sis sp cohe]; % [spatial info (b/s), spatial info (b/sp), sparsity, coherence]
                pdata.(part_now).dwellmap = single(dwellmap);
                pdata.(part_now).speedlift = speedlift;
                if isempty(bdata) || ~any(bdata.partn==sdata_temp.partn) % this part has not been added to bdata yet
                    bdata_temp.dwellmap = { pdata.(part_now).dwellmap };
                end 

%%%%%%%%%%%%%%%% Place fields
                % threshold ratemap
                zmap = ( ratemap - mean(ratemap(:),'omitnan') ) / std(ratemap(:),'omitnan'); % zscore ratemap
                thresh_ratemap = imbinarize(zmap,mapset.zcut); % 2 s.d. threshold
                
                % detect contiguous regions
                datout = regionprops('table',thresh_ratemap,zmap,'Area','Centroid','WeightedCentroid','MajorAxisLength','MinorAxisLength','Orientation','ConvexHull','PixelIdxList','MaxIntensity');
                
                % filter out small or low firing regions
                if ~isempty(datout)
                    datout.Area(:) = datout.Area(:) .* ((mapset.binsize/10)^2); % convert field area to cm2
                    nindx = datout.Area(:) < mapset.arcut | datout.MaxIntensity(:) < mapset.frcut;                        
                    datout(nindx,:) = [];      
                end

                % accumulate data
                sdata_temp.nfields = size(datout.Area,1);
                sdata_temp.field_data = { table2cell(datout) };

%%%%%%%%%%%%%%%% Grid score
                % create autocorrelation
                % I'm not entirely sure who to attribute this to, but autocorrelations are pretty ubiquitous now
                % this function was adapted from xPearson (adapted for speed and to make it work on n-dimensions)
                automap = ndautoCORR(ratemap,ratemap,50);

                % autocorrelation analysis
                % There are a number of different grid score approaches contained here so see the function itself for a description
                % my preferred method is the default, which is from Langston et al. (2010) Development of the Spatial Representation System in the Rat
                % In this method rings are cut from the autocorrelation at different distances from the centre and the standard rotation and correlation
                % is performed on each one. The highest of these grid scores is used, the radius of the ring is the grid spacing. A sine wave
                % is fitted to the values in the ring and this is used to estimate the positions of the fields. The grid orientation is the
                % angle from the centre to the first one of these fields, counter-clockwise. The grid field radius is taken as the radius of the centre field.
                [~,gdata] = get_grid_score(automap,mapset.binsize/10,'method','allen');

                % accumulate data
                sdata_temp.gridmap = { single(automap) };
                sdata_temp.grid_info = single([gdata.grid_score gdata.wavelength gdata.grid_orientation]);
                sdata_temp.grid_field_info = single([gdata.radius gdata.majaxislength gdata.minaxislength gdata.height gdata.width gdata.field_orientation]); 

%%%%%%%%%%%%%%%% HD and Rayleigh vector length
                % HD firing rate maps are now computed by mapHD, like the above it will accept a precomputed
                % dwell time map to speed up computation
                % rayleigh vector values etc are computed on sum normalised HD firing rate maps to remove
                % any possible effects of firing rate (it's unclear if circ_stats care about the values)
                % If you recorded HD using two LEDs this should be used here and throughout the data
                % if not, HD is estimated using the animal's displacement
                % Load a head direction dwell map if one exists, or create it if not
                hd_dwellmap = [];
                if isfield(pdata,part_now)
                    if isfield(pdata.(part_now),'hd_dwellmap')
                        hd_dwellmap = pdata.(part_now).hd_dwellmap;
                    end
                end      
                
                ppoh_map = ppoh;   
                psph_map = psph;
                if mapset.hd_displace
                    ppoh_map = ppod;
                    psph_map = pspd;                    
                end
                [~,hd_dwellmap,hd_spikemap,hd_ratemap,r,mx,mn,sd] = mapHD(hd_dwellmap,ppoh_map,psph_map,mapset);

                % accumulate data
                sdata_temp.hd_ratemap = { single(hd_ratemap) };               
                sdata_temp.hd_info = single([r mx(1) mn sd]); % Rayleigh v, PFD (angle with highest firing), mean angle (norm circ mean), SD of firing (norm circ stdev)
                pdata.(part_now).hd_dwellmap = hd_dwellmap;
                if isempty(bdata) || ~any(bdata.partn==sdata_temp.partn) % this part has not been added to bdata yet
                    bdata_temp.hd_dwellmap = { pdata.(part_now).hd_dwellmap };
                end 
                
                % half maps
                midpoint = median(ppot(:));
                hd_dwellmap_h1 = [];
                hd_dwellmap_h2 = [];                
                if isfield(pdata,part_now)
                    if isfield(pdata.(part_now),'hd_dwellmap_half1')
                        hd_dwellmap_h1 = pdata.(part_now).hd_dwellmap_half1;
                        hd_dwellmap_h2 = pdata.(part_now).hd_dwellmap_half2;                        
                    end
                end      
                [~,hd_dwellmap_h1,~,hd_ratemap_h1,r_h1,mx_h1,mn_h1,sd_h1] = mapHD(hd_dwellmap_h1,ppoh_map(ppot<midpoint),psph_map(pspt<midpoint),mapset); % first half HD map
                [~,hd_dwellmap_h2,~,hd_ratemap_h2,r_h2,mx_h2,mn_h2,sd_h2] = mapHD(hd_dwellmap_h2,ppoh_map(ppot>midpoint),psph_map(pspt>midpoint),mapset); % second half HD map

                % accumulate data
                sdata_temp.hd_ratemap_half = { single(hd_ratemap_h1) single(hd_ratemap_h2) };               
                sdata_temp.hd_info_half = single([r_h1 mx_h1(1) mn_h1 sd_h1 r_h2 mx_h2(1) mn_h2 sd_h2]); % Rayleigh v, PFD (angle with highest firing), mean angle (norm circ mean), SD of firing (norm circ stdev)
                pdata.(part_now).hd_dwellmap_half1 = hd_dwellmap_h1;
                pdata.(part_now).hd_dwellmap_half2 = hd_dwellmap_h2;    
                if isempty(bdata) || ~any(bdata.partn==sdata_temp.partn) % this part has not been added to bdata yet
                    bdata_temp.hd_dwellmap_half = { pdata.(part_now).hd_dwellmap_half1 pdata.(part_now).hd_dwellmap_half2 };
                end 

%%%%%%%%%%%%%%%% Spatial info, grid score, Rayleigh vector shuffles               
                % spatial shuffles
                if config.spatial_shuffles
                    sinfo = savelli_shuffle(ppox.*10,ppoy.*10,ppot,ppoh,pspt,sidx,rmset,100);
                    obs_vals = [sdata_temp.spatial_info(1) mean(table2array(sinfo(:,2)),1,'omitmissing') sdata_temp.hd_info(1)];
                    [z,p] = computeZProbability(obs_vals,table2array(sinfo(:,4:6)));
                    sdata_temp.spatial_info_z = z; % [spatial info (z), grid score (z), Rayleigh vector (z)]
                    sdata_temp.spatial_info_p = p; % [spatial info (p), grid score (p), Rayleigh vector (p)]
                else
                    sdata_temp.spatial_info_z = NaN(1,3); % [spatial info (z), grid score (z), Rayleigh vector (z)]
                    sdata_temp.spatial_info_p = NaN(1,3); % [spatial info (p), grid score (p), Rayleigh vector (p)]
                end                

%%%%%%%%%%%%%%%% AHV analysis
                ahv_window = 200;
                ahv_binsize = 2;
                ahv_mindwell = (1/pos_srate);
                ahv_smoo = 0;
                edg = -ahv_window : ahv_binsize : ahv_window;
                xi = movmean(edg,2,'EndPoints','discard');
                ahv_dwellmap = histcounts(ppoa,edg);
                if ahv_smoo>0
                    ahv_dwellmap = imgaussfilt(ahv_dwellmap,ahv_smoo,'Padding','replicate');
                end
                ahv_dwellmap(ahv_dwellmap<ahv_mindwell) = NaN;
            
                % bin spike ahv into ahv_binsize bins
                ahv_spikemap = histcounts(pspa,edg);
                if ahv_smoo>0
                    ahv_spikemap = imgaussfilt(ahv_spikemap,ahv_smoo,'Padding','replicate');
                end
            
                % ahv ratemap
                ahv_ratemap = ahv_spikemap ./ (ahv_dwellmap .* (1/pos_srate));
            
                % accumulate data
                sdata_temp.ahv_curve = { single(ahv_ratemap(:)) };             
                pdata.ahv_xvalues = single(xi);    
                pdata.(part_now).ahv_dwell = single(ahv_dwellmap);
                if isempty(bdata) || ~any(bdata.partn==sdata_temp.partn) % this part has not been added to bdata yet
                    bdata_temp.ahv_dwellmap = { pdata.(part_now).ahv_dwell };
                end 
            
%%%%%%%%%%%%%%%% Theta phase preference
                if ~isempty(lfp)
                    % bin the theta phase data
                    ai = reshape(deg2rad(-180:5:540),[],1); % doing this means we have bins symmetrical around zero

                    % we must calculate the spike phase here precisely for every spike time, rather than indexing into the position theta phase
                    % this is for reasons of resolution (LFP and spikes are sampled at 48kHz positions are only sampled at 50Hz)
                    % This process is fairly well established:
                    % "To obtain a theta phase angle for each spike, LFPs were first bandpass filtered (fourth-order Chebyshev, r = 0.5, MATLAB 
                    % filter and filtfilt routines; 6–10 Hz) before a Hilbert transform was applied to obtain the instantaneous phase angle."
                    % from van der Meer and Redish (2011) Theta Phase Precession in Rat Ventral Striatum Links Place and Reward Information (https://doi.org/10.1523/JNEUROSCI.4869-10.2011)              
                    pspp = interp1(lfp(:,2),lfp(:,3),pspt,'nearest');
                    spp2 = [pspp(:)-2*pi; pspp(:); pspp(:)+2*pi]; % repeat phase values at +/-360 degrees
                    y = histcounts(spp2,ai); % bin phase values

                    % directional analyses on phase data
                    mx2p = ai(y == max(y)); % preferred angle (location of max frate)

                    % accumulate data
                    sdata_temp.theta_phase = { single(y) }; 
                    sdata_temp.theta_info = single([circ_r(pspp(:)) mx2p(1) circ_mean(pspp)]); % Rayleigh v, preferred phase (bin with max count), mean angle (average phase) 
                else
                    sdata_temp.theta_phase = { NaN }; 
                    sdata_temp.theta_info = NaN(1,3,'single');                   
                end
                
%%%%%%%%%%%%%%%% ISI and autocorrelation analyses
                %% inter-spike interval
                % [~,idata,isis] = getISIhalfwidth(pspt,100);
                [f,xi,e,edg] = compute_hoisa(pspt,'spt2',pspt,'binsize',1,'window',50,'method','isih');

                % accumulate data
                sdata_temp.isi = { f }; 
                sdata_temp.isi_info = single([ NaN ]); 
                pdata.isi_xvalues = single(xi);

                %% burst index
                % from Mizuseki et al. (2012) Activity Dynamics and Behavioral Correlates of CA3 and CA1 Hippocampal Pyramidal Neurons 
                % "The spike-burst index was defined as the fraction of spikes with <6 ms ISIs (Harris et al., 2001)"
                % https://doi.org/10.1002/hipo.22002  
                isis = diff(pspt) .* 1e3; % ISIs in ms
                bindx = [isis; NaN] < 6 | [NaN; isis] < 6; % bindx is an index of all spikes sharing an isi less than 6ms
                burst_index = sum(bindx) / numel(pspt);

                % accumulate data
                sdata_temp.burst_index = single(burst_index); 

                %% theta   
                [f,xi,e,edg] = compute_hoisa(pspt,'spt2',pspt,'binsize',10,'window',500,'method','hoisa');

                % accumulate data
                sdata_temp.autocorr_theta = { single(f(:)) }; 
                % sdata_temp.autocorr_500_info = single([NaN NaN NaN]); % theta index, theta frequency, fit rsquare
                pdata.autocorr_theta_xvalues = single(xi);
                pdata.autocorr_theta_evalues = single(edg);

                %% Refractory period analyses           
                % refractory period violation analyses
                % from Navratilova et al. (2016) Grids from bands, or bands from grids? An examination of the effects of single unit contamination on grid cell firing fields.
                % although they took these analyses from pre-existing papers
                % "There are few reliable methods for estimating the contamination of a unit isolated from tetrode recordings (see DISCUSSION). The only method that does not rely 
                % on the same measures that are used for spike sorting is to check the number of spikes that occur within the refractory period of another spike..."
                % https://doi.org/10.1152/jn.00699.2015
                % [nrpv,prpv,fp1,fp2,censored,t25,c25] = getRPVcount(pspt,isis,part_config.part_duration(pp)); % we extract the number of spikes in the refractory period, what proportion of the total this is etc                    
                [f,xi,e,edg] = compute_hoisa(pspt,'spt2',pspt,'binsize',1,'window',50,'method','hoisa');
                nrpv = sum(isis < 2);
                prpv = nrpv ./ (numel(pspt)-1);

                % accumulate data
                sdata_temp.autocorr_refrac = { single(f(:)) }; 
                sdata_temp.autocorr_refrac_info = single([nrpv prpv]);     
                pdata.autocorr_refrac_xvalues = single(xi);  
                pdata.autocorr_refrac_evalues = single(edg);

%%%%%%%%%%%%%%%% Speed cell analysis
                % This analysis is taken from Kropff et al (mEC speed cell paper)
                % It is unlikely to be of great use to HPC people, although many place cells are speed modulated
                % Related to this, you can also compute the relationship between theta power and running speed
                [svals,stime,sscore,sslope,sintcpt,crve,~] = getSPEEDmod(ppot,ppov,pspt);

                % accumulate data
                sdata_temp.speed_slope = { single(interp1(svals,crve,0:1:50,'linear')) };             
                sdata_temp.speed_info = single([sscore sslope sintcpt]);             
                pdata.(part_now).speed_dwell_time = single(interp1(svals,stime,0:1:50,'linear'));
                if isempty(bdata) || ~any(bdata.partn==sdata_temp.partn) % this part has not been added to bdata yet
                    bdata_temp.speed_dwell_time = { pdata.(part_now).speed_dwell_time };
                end 
                
%%%%%%%%%%%%%%%% Accumulate data
                sdata = [sdata; sdata_temp];
                bdata = [bdata; bdata_temp];
                
                % Create summary figure for this cluster in this part
                fname = fullfile(pwd,pdata.outname,'part_figures',[sdata_temp.uci{1} '_' part_now '.png']);
                if exist(fname,'file') && config.skipfigs % if the figure exists and we don't want to overwrite it
                    continue % don't make the figure
                elseif isempty(pspt) || ~config.save_figs % if there were no spikes in this part or we don't want figures at all
                    continue % don't make the figure  
                else
                    % putvar(pdata,sdata_temp,pp,wav_now_clus,quals,fets,fast_figures); return;
                    klustfig_part(pdata,sdata_temp,pp,wav_now_clus,quals,fets,config);
                end
            end
        end
    end
    
%%%%%%%%%%%%%%%% Save the data and finish up
    sdata = addprop(sdata,{'pdata'},{'table'});
    sdata = addprop(sdata,{'bdata'},{'table'});    
    sdata.Properties.CustomProperties.pdata = pdata;
    sdata.Properties.CustomProperties.bdata = bdata;   
    fname = fullfile(pwd,pdata.outname,'sdata.mat');
    save(fname,'sdata'); % save session data
    analysis_log({'klustest'},1,'version',{'v20.0.4'});
    
    % finish up
    toc1 = toc/60;
    disp(sprintf('klustest has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',fullfile(pwd,pdata.outname),' &'');">','klustest folder','</a>'])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');

