function klustest(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
%klustest  function for analysing clusters after cluster cutting
%    This function utilises many minor functions to both analyse and plot cluster data generated in Tint. If requested, 
%    it will output a figure for each cluster, the cluster space of each tetrode, the cross-correlations of every cluster 
%    on a tetrode and a session data structure (sdata.mat) It will also generate an mtint file (mtint.mat) containing all 
%    the tetrode and cluster info.
%
% USAGE:
%           klustest() process files with default settings
%                 klustest will look for .cut files named [cname] (default is kwiktint), it will open the first one and extract 
%                 details about which sessions were used to make that .cut file. It will then open all of the required files to
%                 build an mtint table which contains a summary of all the spikes etc for a session. It will then process every
%                 cluster on every possible tetrode and generate an output figure for each.
%
%           klustest(Name,Value,...) process with Name-Value pairs used to control aspects of the function output
%
%           Parameters include:
%
%           'tetrodes'              -   (default = 1:16) Vector, the tetrodes to run on (i.e. [1 2 3 4 5 6] would run on tetrodes 1 to 6)
%
%           'clusters'              -   (default = 0) Vector, the clusters to run on, set this to 0 to run on all clusters
%
%           'rname'                 -   (default = first 3 characters of parent directory) String, the rat name/number; this will be used in the sdata structure so its important to give this
%
%           'cname'                 -   (default = 'kwiktint') String, function will automatically look for klustakwiked files named this
%
%           'part_names'            -   Cell array of strings, the names you want the outputs to be saved as - see function text for more info
%                                       default values are whatever values are specified in klustest
%
%           'part_methods'          -   Numeric vector, see function text for more info
%                                       default values are whatever values are specified in klustest
%
%           'part_intervals'        -   Cell array of numeric vectors, see function text for more info
%                                       default values are whatever values are specified in klustest
%
%           'pixel_ratio_override'  -   (default = []) Numeric vector, set to empty [] to ignore, otherwise this must be n recording sessions long
%
%           'interval_keys'         -   (default = {'s','e'}) Cell array, what keypresses were used to delineate trial starts and ends {start,end} these can be integers or characters - see function text for more info
%
% EXAMPLES:
%
%           % run function using default values
%           klustest()
%
%           % run function using default values, but only on tetrodes 1 and 5
%           klustest('tetrodes',[1 5])
%
%           % run function using default values, all specified
%           klustest('tetrodes',1:16,'clusters',0,'rname','RAT','cname','kwiktint')
%
% See also: kwiktint getTRODESkk

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

%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2016 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INPUT ARGUMENTS CHECK
%% Prepare default settings
    def_tetrodes                = 1:16;        
    def_clusters                = 0;     
    def_cname                   = 'kwiktint';   
    % I can extract the rat name from the file path because I know my folder structure, if this is not the case for you, change def_rname to what you prefer
    % the next two lines make rname equal the first 3 characters of the directory containing the current directory
    pname = pwd; sindx=strfind(pname,'\'); 
    def_rname                   = pname(sindx(end-1)+1:sindx(end-1)+3);    
    
    % Specify partition configuration settings
    % If these are passed as name,value pair arguments those values will be taken instead
    def_part_names              = {'S1'}; % the names you want the outputs to be saved as, these cannot start with a number: i.e. 'session1'
    % method of partition, corresponding to each of the names above: 1 = combine everything, 2 = take a recording session, 3 = use digital inputs
    def_part_methods            = [2]; % i.e. if part_methods=[2 2 2 1], we will have 4 parts, the first 3 correspond to some combination of recording sessions (there can be multiple ones) and the last one will include all data
    % cell array of vectors indicating which intervals to include in each partition: if method = 1 this does nothing, if method = 2 this should specify which recording sessions to include, if method = 3 this should specify which digital input pairs to include (inf = use all)
    def_part_intervals          = {1}; % i.e. if part_methods=[2 2 2], then if part_intervals={1 2 3}, rec 1 will go in part 1, rec 2 in part 2 and rec 3 in part 3    OR     if part_methods=[2 2 2 2 2], then if part_intervals={1 [2 4 5] 3}, rec 1 will go in part 1, rec 2,4 and 5 in part 2 and rec 3 in part 3
    def_pixel_ratio_override    = []; % set to empty [] to ignore, otherwise this must be n recording sessions long
    def_interval_keys           = {{'s','e'}}; % what keypresses were used to delineate trial starts and ends {start,end} these can be integers or characters i.e. {1,2} or {'s','e'}, {'1','2'} is the same as {1,2}
            
%% Parse inputs
    p = inputParser;
    addParameter(p,'tetrodes',def_tetrodes,@(x) ~isempty(x) && ~all(isnan(x(:))) && isnumeric(x));  
    addParameter(p,'clusters',def_clusters,@(x) ~isempty(x) && ~all(isnan(x(:))) && isnumeric(x)); 
    addParameter(p,'cname',def_cname,@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    addParameter(p,'rname',def_rname,@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    addParameter(p,'part_names',def_part_names,@(x) ~isempty(x) && iscell(x)); 
    addParameter(p,'part_methods',def_part_methods,@(x) ~isempty(x) && isnumeric(x)); 
    addParameter(p,'part_intervals',def_part_intervals,@(x) ~isempty(x) && iscell(x)); 
    addParameter(p,'pixel_ratio_override',def_pixel_ratio_override,@(x) ~isempty(x) && isnumeric(x)); 
    addParameter(p,'interval_keys',def_interval_keys,@(x) ~isempty(x) && iscell(x) && numel(x)==2); 
    parse(p,varargin{:});

%% Retrieve parameters 
    config = p.Results;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INITIAL SETTINGS / INPUTS
    % overrides - general settings for overriding normal klustest functionality
    config.pconfig_override = 0; % set to 1 if you want to ignore and overwrite an existing part_config, set to 2 to run with current part_config settings without overwriting anything
    config.maintain_mtint   = 0; % DEBUGGING ONLY set to 1 to save/load mtint in the base workspace, this saves time when running the function mutliple times (for instance in debugging) but should otherwise be set to 0
    config.mtint_override   = 0; % set this to 1 to force the production of a new mtint file (do this if you reclustered in tint for instance)
    config.trial_override   = 0; % set this to 1 to force the production of a new .key file (do this if you are not happy with the current trials for instance)
    config.fig_key_off      = 1; % set this to 1 to suppress display of trials/digital keypresses in an interactive UI
    
    % map settings - settings used when generating 2D spatial firing rate maps
    config.rmethod          = 'gaussian'; % (default 'nearest') the mapping approach to use, either 'nearest','gaussian','adaptive','KDE'
    config.map_padd         = 2; % (default 2) the number of bins to pad spatial maps with
    config.bin_size         = 2; % (default 2) bin size in cm for calculating the rate map (make sure data is in cm or pm/sm values are given)
    config.map_sigma        = 1.5; % (default 1.5) used by nearest and KDE method, sigma of gaussian to use when smoothing traditional dwell and spike maps, or used as bandwidth of KDE
    config.smethod          = 1; % (default 1) used by nearest method, when set to 1 the dwell and spike maps are smoothed before dividing, 2 means the smoothing is applied after dividing, 3 means no smoothing
    config.min_dwell        = 0.01; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    config.g_sigma          = 5; % (default 10) only used by gaussian method - sigma of gaussian used to weight position point and spike point distances, bigger means smoother maps
    config.min_dist         = 5; % (default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    config.max_dist         = 20; % (default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    config.srate            = 50; % (default 50) sampling rate of data in Hz, used to calculate time

    % field settings - settings to use when detecting place fields
    config.frcut            = 0.5; % relative minimum firing rate (% of ratemap max) to be considered a field
    config.arcut            = 9; % minimum number of pixels to be considered a field, typically people use 9 contiguous pixels
    config.minfr            = 1; % (Hz) absolute minimum cutoff firing rate to be considered a field

    % spike plot settings
    config.over_smooth      = 13; % number of position data points over which to smooth instantaneous firing rate when calculating overdispersion
    config.time_bins        = 2; % (default 2s) time window over which to compute the spike vs time plot

    % HD settings - settings to use when generating head direction maps
    config.hd_type          = 'histogram'; % (default 'histogram') enter 'density' for a circular kolmogorov smirnov density estimate plot, enter 'histogram' for the traditional histogram polar plot
    config.hd_bins          = 64; % (default 64) the number of bins to use when computing HD plot
    config.hd_sigma         = 0.04; % (default 2) the standard deviation of the gaussian kernel used in HD circular density estimate
    config.hd_boxcar        = 3; % (defualt 3) the number of bins over which to compute the HD histogram boxcar average

    % figure settings
    run_figPART            = 1; % [ figPART ] set to 1 for figures showing the activity of a cluster in each part
    run_figCLUS            = 0; % [ figCLUS ] set to 1 for figures showing the activity of a cluster in all parts together
    run_figCROSS           = 0; % [ figCROSS ] set to 1 for figures showing spike cross-correlation of all cells on a given tetrode   
    fig_vis                = 'off'; % set to 'on' to see figures as they are generated and saved, set to 'off' to suppress this (faster)
    
    % optimization settings
    config.wave_asis        = 1; % set to 1 for waveform analysis and plots
    config.temp_asis        = 1; % set to 1 for spike time analysis and plots
    config.fild_asis        = 1; % set to 1 for place field analysis and plots
    config.grid_asis        = 1; % set to 1 for grid cell analysis and plots
    config.head_asis        = 1; % set to 1 for HD cell analysis and plots
    config.ovrd_asis        = 1; % set to 1 for overdispersion analysis and plots
    config.sisi_asis        = 1; % set to 1 for inter-spike interval analysis and plots
    config.autc_asis        = 1; % set to 1 for spike autocorrelation analysis and plots
    config.spph_asis        = 1; % set to 1 for spike theta phase analysis and plots
    config.sped_asis        = 1; % set to 1 for speed modulation analysis and plots
    config.cell_asis        = 1; % set to 1 for cell type analysis and plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Start analysis and get part_config
    % get the function's name (people may rename it from klustest)
    stk = dbstack;
    function_name = stk.name;

    % starting messages
    tic;
    disp('----------------------------------------------------------------------------');
    time_now = datestr(now,'yyyy-mm-dd-HH-MM-SS');
    disp(sprintf('Running %s at %s...',function_name,time_now))
    if config.pconfig_override; disp(sprintf('\tWARNING: will override part_config...'));                 end
    if config.maintain_mtint;   disp(sprintf('\tWARNING: will use/maintain mtint in base workspace...')); end
    if config.mtint_override;   disp(sprintf('\tWARNING: will override mtint...'));                       end
    disp(sprintf('\t...%s spatial maps with %.1fcm bins (%.1f padding)',config.rmethod,config.bin_size,config.map_padd));
    disp(sprintf('\t...%.1f%% cutoff for fields, %.1f pixels minimum, %.1fHz minimum',config.frcut.*100,config.arcut,config.minfr));
    disp(sprintf('\t...%s HD maps with %d bins',config.hd_type,config.hd_bins));

    % The part_config table contains basic information about how/where/if to separate the data into different
    % partitions or parts. For instance if we record an open field, then some sort of maze, then the open field
    % again, we may want to divide the data into those 3 parts for separate analysis, or we might want to just 
    % lump everything together if we care about some intrinsic property of the cells
    % We want to save this table so that in the future we can just run the function with the same settings
    [~,~,~] = mkdir([pwd '\klustest\' config.cname]); % create the folder which will hold a lot of klustest data
    part_num = length(config.part_names);
    part_config = table(repmat({config.rname},part_num,1),config.part_names',config.part_methods',config.part_intervals',repmat({config.cname},part_num,1));
    part_config.Properties.VariableNames = {'Rat','Part','Method','Intervals','Outname'};
    
    % load or save part_config
    % as mentioned above, the part_config file contains the basic information required to process a dataset, it really
    % is just a permanent copy of the settings outlined manually at the top of klustest. For more information, like
    % the start and end times of each part, or the files used for each one, see the part_config saved in pdata instead
    % partIO here will try to save the part_config in a file named pconfig_name, if this file already exists it will
    % backup the contents using the current date/time and then save the part_config alongside it non destructively
    % Otherwise it will make a new file and save it there. FOr this reason, if you load pconfig_name you will actually find 
    % a structure named 'part_data', the field named 'part_config' is the last part_config setting used, the other fields are
    % older part_configs, where the numbers in the field name reflect the date/time they were overwritten
    pconfig_name = ['klustest\' config.cname '\' config.cname '_part_config.mat'];    
    part_config = partIO(pconfig_name,part_config,config.pconfig_override);
    disp(sprintf('\t...rat name: %s',config.rname))      

    % start to prepare the pdata (part or session data) and sdata (cell data) tables
    % To save space and make things easier for the user I have divided the data into two files
    % the pdata structure contains data which applies to the whole session, like positions or dwell maps
    % the sdata table contains data for each cluster/part (one per row). The idea being that tables are easier to concatenate
    % to build a full data set
    pdata = struct; % will hold part data
    pdata.part_names = part_config.Part';
    pdata.combined_name = part_config.Outname{1,1};
    pdata.rat = part_config.Rat{1,1};
    pdata.analysed = time_now;
    pdata.directory = pwd;
    pdata.interval_keys = config.interval_keys;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Get session names and tetrodes
    % When sessions are cluster cut together (which you should do with multiple sessions recording the same cells)
    % tint saves the filenames of the individual sessions in the .cut file it makes after cluster cutting with
    % klustakwik. I found extracting this data to be the easiest way to get this information, and it ensures the
    % function will only ever be used to analyse sessions cluster cut in one file
    % As it happens, we can also check to see which tetrodes have a cut file and from this work out which tetrodes
    % we can and can't analyse
    disp(sprintf('Assessing data...'))
    [tets,snames,cutname] = getTRODES(config.cname); 
    
    disp(sprintf('\t...read %s',cutname));
    nsess = numel(snames);
    disp(sprintf('\t...working on %d sessions: %s',nsess,strjoin(snames,', ')));    
    if ~isempty(config.tetrodes) % if specific tetrodes were requested
        if ~isempty(setdiff(config.tetrodes,tets)) % if there were some tetrodes requested that don't have files
            disp(sprintf('\t...WARNING: skipping tetrodes %s, cut file not found',mat2str(setdiff(config.tetrodes,tets))))            
        end
        tetrodes = intersect(config.tetrodes,tets);
    else % if the user is running in automatic detection of tetrodes
        tetrodes = tets;
    end
    disp(sprintf('\t...tetrodes: %s accounted for',mat2str((tetrodes(:))')))

    % accumulate
    pdata.tetrodes = tetrodes; % list of tetrodes analysed
    pdata.session_names = snames; % names of the recording sessions (i.e. files) used, cell array
    pdata.sessions = nsess; % the number of recording sessions used
    pdata.part_config = part_config; % a copy of the part_config, this will be extended to include more data that the saved part_config though, table format
    pdata.config = config; % a copy of the configuration struct, so ratemaps etc can be remade with identical settings if need be

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Load or generate mtint
    % The first step of the data analysis is just to load all of the Axona file data into a matlab friendly format
    % this is largely unprocessed, raw data. This is stored in a structure called mtint.
    % We create this and will use it a little bit, but mostly it just serves as a foundation for klustest. Our 
    % data and stats will be stored in the pdata and sdata structures. Klustest shouldn't really change the contents
    % of this structure as it serves as a non-destructive reference
    % Because the mtint can be quite large and bulky, we have some time saving options - if one has already been created 
    % we can just load it instead of making a new one. If you re-cluster cut the data you should override any existing
    % mtint though (or delete it manually) otherwise the clusters will not match those in Tint. If we are debugging the function
    % we can just opt to use an mtint stored in the base workspace (so no loading necessary)
    disp(sprintf('Fetching DACQ data (mtint)...'));
    mname = ['klustest\' config.cname '\' config.cname '_mtint.mat']; % the filename of the mtint (one that exists or one we will make)
    
    if any(strcmp(evalin('base','who'),'mtint')) && ~config.mtint_override && config.maintain_mtint % if we want to use an mtint currently held in memory
        mtint = evalin('base','mtint');        
        m = whos('mtint');     
        disp(sprintf('\t...using mtint held in memory (%.1fMb)',m.bytes./1e+6));
    elseif ~exist(mname,'file') || config.mtint_override
        disp(sprintf('\t...building mtint'));    
        mtint = getDACQDATA(config,snames,tetrodes); % my replacement for readalldacqdata
        disp(sprintf('\t...post-processing mtint'));

        % deal with manual override of pixel ratio if necessary
        if numel(config.pixel_ratio_override) ~= numel(mtint.pos.header)
            disp(sprintf('\tWARNING: number of pixel ratios in part_config (%d) does not equal number of recordings (%d)... ignoring manual values',numel(config.pixel_ratio_override),length({mtint.pos.header.pixels_per_metre})));
        else
            for pr = 1:numel([mtint.pos.header.pixels_per_metre])
                mtint.pos.header(pr).pixels_per_metre = config.pixel_ratio_override(pr);
            end
        end    

        % post process position data etc
        % this function adds running speed, head direction data and also smoothes, interpolates the data
        % removes jumps in the data which are too quick to be natural, converts the position data to cm
        % and a bunch of other stuff
        mtint = postprocessDACQDATA(mtint); % my replacement for postprocess_DACQ_data

        % save mtint
        info = whos('mtint');
        siz = info.bytes / 1000000;
        disp(sprintf('\t...saving mtint (%.1fMb)',siz)) % I have tried my best to keep the mtint >50Mb for most recordings
        save(mname,'mtint','-v7.3');
    else
        m = dir(mname);
        disp(sprintf('\t...loading mtint (%.1fMb)',m.bytes./1e+6));
        load(mname,'mtint');
    end

    % basic session info
    duration = mtint.pos.total_duration;
    disp(sprintf('\t...total session time: %ds',duration))
    disp(sprintf('\t...cut file is made up of %d recordings',numel(mtint.header)))

    if config.maintain_mtint
        assignin('base','mtint',mtint); % leave mtint in base workspace
    end 
    disp(sprintf('\t...done'));

    % accumulate
    pdata.date = mtint.header_date; % date of recording
    pdata.duration = duration; % total duration of the recording (s)
    pdata.recording_times = mtint.pos.trial_duration; % the length in s of each recording
    pdata.session_details = table; % table to hold extra session info in one place
    pdata.session_details.session_names = pdata.session_names; % cell array of session names
    pdata.session_details.session_lengths = mtint.pos.trial_duration(:); % length of each recording (s)
    pdata.session_details.session_ends = cumsum(mtint.pos.trial_duration(:)); % the time (in concatenated klustest time) where this session ends in the data (s)
    pdata.session_details.session_starts = pdata.session_details.session_ends - pdata.session_details.session_lengths; % the time (in concatenated klustest time) where this session starts in the data (s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Add interval/part time delineations to part_config (now stored in pdata)
    % We want to divide a long recording session into partitions or parts, these might be whole
    % recording sessions or they might be based on digital keypresses or you might just want the
    % whole thing in one long session. In any case partSESS2 takes the keypresses and the desired
    % method defining each part and gives the start + end time of each so we can separate the data later
    % If you are using keypresses (i.e. intervals) partSESS2 will also run figKEYS which will plot every
    % trial and let you correct keypresses if necessary
    disp(sprintf('Preparing partitions...'))
    pdata.config.trial_override = config.trial_override;
    pdata.config.fig_key_off = config.fig_key_off;    
    [pdata,nparts] = partSESS2(mtint,pdata); % process the part_config to get the start and end time of each partition

    % basic session info
    disp(sprintf('\t...%d parts',nparts))
    disp(sprintf('\t...part methods: %s',mat2str(pdata.part_config.Method)))
    disp(sprintf('\t...part names: %s',strjoin(pdata.part_config.Part,', ')))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Prepare important data
    disp(sprintf('Extracting initial data...'))
    % get the position data (dacqUSB data) for the whole session
    % position data are already converted to cm and saved in the mtint structure so we will just load it here
    % most stuff in the mtint is saved as single so we also need to convert that out (some operations require double)
    % also saved in the mtint are the head direction data (or displacement if only 1 LED was used) and pixels per metre
    % converted into column format to match the position data
    pox = double(mtint.pos.xy_cm(:,1)); % extract just the x coordinates
    poy = double(-mtint.pos.xy_cm(:,2)); % extract just the y coordinates, reflect these as dacqUSB tracks from top left corner of picture
    pot = double(mtint.pos.ts(:)); % extract the time stamp of each position value
    pov = double(mtint.pos.speed(:)); % the running speed throughout the session
    poh = double(mtint.pos.dir(:,1)); % HD data
    ppm = double(mtint.pos.pixels_per_metre(:)); 

    % get the position data sampling rate (should be 50hz) or 0.05s
    srate = mtint.pos.header(1).sample_rate_num(1,1); % sampling rate
    sinterval = 1 / srate; % sampling interval

    % display results
    disp(sprintf('\t...positions read: %d',numel(pox(:))));
    disp(sprintf('\t...tracking LEDs: %d',size(mtint.pos.led_pos,2)));
    disp(sprintf('\t...median pixel ratio: %dppm',nanmedian(ppm)));
    disp(sprintf('\t...sample rate: %dHz (%.2fs)',srate,sinterval));

    % get the LFP data
    % as explained below in the LFP section, klustest loads the first available LFP data saved in mtint, this corresponds to
    % LFP channel 1 in dacqUSB (i.e. the primary LFP channel displayed during recording)
    % this can be changed here if necessary, but if we try to load a non-existant field of mtint the function will crash
    Fs = mtint.lfp(1).Fs(1,1);
    lfp = double(mtint.lfp(1).lfp(:,1));
    lfpt = (0:length(lfp)-1)'/Fs; % make a vector for time    

    % filter LFP to get the theta band
    [b,a] = butter(4,[6 12]/(Fs/2)); % Generate 4th order butterworth filter coefficients for theta band [6 12] Hz
    lftheta = filtfilt(b,a,lfp); % Apply filter to data using zero-phase filtering  
    
    % accumulate data
    pdata.pox = pox;
    pdata.poy = poy;
    pdata.pot = pot; 
    pdata.pov = pov;
    pdata.poh = poh;    
    pdata.ppm = ppm; 
    pdata.sampling_rate = srate;
    pdata.sampling_interval = sinterval;
    pdata.leds = size(mtint.pos.led_pos,2);  

    % display results
    disp(sprintf('\t...LFP samples read: %d',numel(lfp)));
    disp(sprintf('\t...LFP sample rate: %dHz',Fs));    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################################################################################## %% Run through tetrodes and clusters
    disp(sprintf('Analysing tetrode/cluster data...'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% For every available tetrode
    % The next loop focuses on analysing a single tetrode, ther eis not a lot of info we really need about tetrodes
    % because the meat of the data is really per cluster or part
    % This is especially true now that cluster quality has been moved to gatDACQDATA, so that it doesn't have to be
    % computed on every run of klustest. However, we will load the waveforms for this tetrode before continuing to the clusters
    sdata = table; % will hold all the cluster data
    bdata = table; % will hold all the behaviour data (most of it duplicated from pdata though)
    for ee = 1:length(tetrodes) 
        tet = tetrodes(ee); % tet = the current tetrode
        disp(sprintf('\tLoading tetrode %d...',tet));

        % get a vector of clusters we want to analyse
        if isempty(config.clusters) || config.clusters==0
            clus = unique(mtint.tetrode(tet).cut); % find the clusters logged for this tetrode
        else
            clus = config.clusters;
        end
        pdata.clusters{tetrodes(ee)} = clus;

        % check to see if there are any clusters
        clus_count = numel(clus);   
        disp(sprintf('\t\t...%d spikes detected',mtint.tetrode(tet).nspike_cut));
        disp(sprintf('\t\t\t...%d data clusters detected',sum(clus~=0))); 
        if ~clus_count || ~any(clus) % if there are no clusters, or if there is only a noise cluster
            continue % skip analysis
        end
        disp(sprintf('\t\t\t\t...starting analysis'));

        % get the channel waveforms for this tetrode    
        if config.wave_asis
            disp(sprintf('\t\t\t\t...getting waveforms'));                   
            nspikes = mtint.tetrode(tet).nspike_cut;
            waves = zeros(nspikes,50,4,'single');
            rnow = 1;            
            for ssn = 1:length(pdata.session_names)
                fnamen = pdata.session_names{ssn};
                [~,c1,c2,c3,c4,~] = get_SPIKESV([fnamen,'.',num2str(tet)]);               
                waves(rnow:rnow+length(c1(:,1))-1,:,1) = c1;
                waves(rnow:rnow+length(c2(:,1))-1,:,2) = c2;
                waves(rnow:rnow+length(c3(:,1))-1,:,3) = c3;
                waves(rnow:rnow+length(c4(:,1))-1,:,4) = c4;                
                rnow = rnow+length(c4(:,1));
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% For every detected cluster    
        % The next loop focuses on analysing a single cluster from this tetrode, almost immediately we will skip it if
        % the cluster is noise (i.e. cluster 0), in some cases it might be interesting to have that information, but as we
        % dump all our noise spikes in cluster 0 this may be extraordinarily large, so it makes more sense to skip it for speed

        % get all the spike times for this tetrode
        spiketime = double(mtint.tetrode(tet).ts);      
        disp(sprintf('\t\t\t\t...analysing clusters'));                   
    
        loopout = looper(length(clus));
        for cc = 1:length(clus) % for every cluster   
            sdatac = table; % to hold this cluster's data
            clu = clus(cc); % clu = the current cluster
            uci = [config.rname '.' pdata.date '.' num2str(tet) '.' num2str(clu)]; % unique cell identifier - a string with [rat name, session date, electrode, cluster];

            % if this is the noise cluster don't continue any further
            if ~clu      
                loopout = looper(loopout);
                continue
            end
            
            % clu_indx is a logical vector, length = N of spikes, true if spike is in cluster clu, false if not
            % mtint.tetrode.cut contains a vector specifying which spikes are in which cluster, each row corresponds to a tetrode
            % so to get the data for one tetrode we need to do: mtint.tetrode(tetrode # we want).cut
            clu_identity = mtint.tetrode(tet).cut; % clu_assign is a vector of numbers, one for each spike, each number corresponds to a cluster                        
            clu_indx = clu_identity==clu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% For each part  
            % The next loop focuses on dividing the data into its different parts and analysing the cell within each of these
            % The very first section deals with extracting the data relevant to this part. I've tried to keep this process
            % logical so as to increase speed, the downside is that if there is a large amount of data this approach may fail
            % due to a lack of memory. However, as the matrices are logical (i.e. 1 byte per value) I imagine this would require
            % an enormous recording session with many intervals
            part_names = pdata.part_config.Part;
            for pp = 1:nparts % for every partition         
                sdatap = table; % to hold cluster data for this part  
                bdatap = table; % to hold behaviour data for this part
                part_now = part_names{pp}; % the name of the current part
                if ~isfield(pdata,part_now)
                    continue
                end
                part_times = pdata.(part_now).times; % the time pairs (intervals) corresponding to this part
        
                % find what data falls into the intervals associated with this part
                % start with positions (load existing ones if possible)
                if isfield(pdata.(part_now),'pox')
                    ppox = double(pdata.(part_now).pox); % pos x
                    ppoy = double(pdata.(part_now).poy); % pos y
                    ppot = double(pdata.(part_now).pot); % pos time
                    ppov = double(pdata.(part_now).pov); % pos running speed
                    ppoh = double(pdata.(part_now).poh); % pos HD      
                    part_duration = pdata.(part_now).duration;                     
                else
                    pindax = logical(sum(pot' > part_times(:,1) & pot' < part_times(:,2),1));
                    ppox = pox(pindax); % pos x
                    ppoy = poy(pindax); % pos y
                    ppot = pot(pindax); % pos time
                    ppov = pov(pindax); % pos running speed
                    ppoh = poh(pindax); % pos HD   
                    part_duration = sum(pindax(:))*pdata.sampling_interval; 
                    
                    % accumulate
                    pdata.(part_now).pox = single(ppox);
                    pdata.(part_now).poy = single(ppoy);
                    pdata.(part_now).pot = single(ppot); 
                    pdata.(part_now).poh = single(ppoh); 
                    pdata.(part_now).pov = single(ppov);
                    pdata.(part_now).duration = part_duration;
                end

                % spike data next
                sindax = logical(sum(spiketime' > part_times(:,1) & spiketime' < part_times(:,2) & clu_indx',1)); % spikes
                nindax = logical(sum(spiketime' > part_times(:,1) & spiketime' < part_times(:,2) & ~clu_identity',1)); % noise 
                pspt = spiketime(sindax);
                sidx = knnsearch(pot,pspt);
                sidx2 = knnsearch(ppot,pspt);                
                pspx = pox(sidx); % spike x
                pspy = poy(sidx); % spike y
                psph = poh(sidx); % direction
                pspv = pov(sidx); % running speed                

                % accumulate, these variables don't need to be preallocated as they should always be filled
                sdatap.rat = {pdata.rat};
                sdatap.date = str2double(pdata.date);
                sdatap.directory = {pwd};
                sdatap.partn = pp;
                sdatap.part = {part_now};                
                sdatap.uci = {uci};
                sdatap.tet = tet;
                sdatap.clu = clu;
                sdatap.spike_index = {uint32(sidx)}; % this index can be used to get the spike values from the position data values i.e. pdata.(part_now).pox(sdata.spike_index) would be the spike X positions
                sdatap.part_spike_index = {uint32(sidx2)}; % this index can be used to get the spike values from the position data values i.e. pdata.(part_now).pox(sdata.spike_index) would be the spike X positions                
                sdatap.spike_times = {double(pspt)}; % the actual spike times
                sdatap.isod = mtint.clu_quality(tet).isolation_distances(clu);
                sdatap.lratio = mtint.clu_quality(tet).lratios(clu);
                sdatap.spikes = numel(pspx);
                sdatap.duration = part_duration;
                sdatap.frate = numel(pspx) / part_duration;
                minspikes = 1; % the minimum number of spikes a cluster has to have before we do the following analyses
                pdata.minspikes = minspikes;
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Waveforms
                % all waveforms for this electrode are loaded above by getSPIKESV, but we want only the waveforms
                % corresponding to this cluster and part, which is what this section is for
                % The waveforms are sadly not saved anywhere because they would take up a lot of space (minimum 100mb)
                % so instead we have to just load them from the spike files every time we run klustest
                % Instead we save the mean +/- SD for each channel as a representation of the cluster
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'waveform_mean',cell(1,4),'waveform_stdv',cell(1,4),'waveform_rms',NaN(1,4),'waveform_snr',NaN(1,4),'waveform_max',NaN(1,4),'waveform_min',NaN(1,4),'waveform_maxt',NaN(1,4),'waveform_mint',NaN(1,4),'waveform_width',NaN(1,4),'waveform_params',NaN(1,2),'waveform_snrs',NaN(1,4));

                % get the waveforms for this cluster   
                if config.wave_asis && numel(pspx)>minspikes
                    waves2 = waves(sindax,:,:);
                    for w = 1:4 % for every recording channel
                        wav = double(waves2(:,:,w));
                        ch = nanmean(wav,1);
                        chs = nanstd(wav,[],1);
                        chrms = nanmean(rms(wav,1));
                        nsrms = nanmean(rms(waves(nindax,:,w),1));
                        [maxval,maxt] = max(ch);
                        [postminval,postmint] = min(ch(maxt:end));
                        postmint = postmint + maxt - 1;
                        width = postmint - maxt;
                        width = width * (1000/50); 

                        sdatap.waveform_mean(1,w) = {ch};
                        sdatap.waveform_stdv(1,w) = {chs};
                        sdatap.waveform_rms(1,w) = chrms;  
                        sdatap.waveform_snr(1,w) = chrms/nsrms; % https://en.wikipedia.org/wiki/Signal-to-noise_ratio                       
                        sdatap.waveform_max(1,w) = maxval;
                        sdatap.waveform_min(1,w) = postminval;   
                        sdatap.waveform_maxt(1,w) = maxt;
                        sdatap.waveform_mint(1,w) = postmint;                         
                        sdatap.waveform_width(1,w) = width;
                        
                        % signal to noise Liu et al. (2014) Quality Metrics of Spike Sorting Using Neighborhood Components Analysis
                        % https://dx.doi.org/10.2174%2F1874120701408010060
                        snrs = 1/size(wav,1) .* nansum((nanmax(wav,[],2)-nanmin(wav,[],2)) ./ (2.*nanstd(wav-nanmean(wav,1),[],2)));
                        sdatap.waveform_snrs(1,w) = snrs;
                    end

                    % get the cell's width of waveform using the waveform with the highest mean amplitude
                    [amp,idx] = max(sdatap.waveform_max);
                    wow = sdatap.waveform_width(1,idx);
                    sdatap.waveform_params(1,:) = [amp wow]; % accumulate data
                    clear waves2 sindax pindax % clear all the waveform data, we don't need it again and it is quite large
                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Ratemap, dwellmap, spatial analyses, grid analyses, overdispersion
                % place field, gridness and overdispersion analyses all require a firing rate map
                % so they are combined here in the same loop where the firing rate map is generated
                % firing rate maps are generated by mapDATA which can utilise a number of different methods
                % it will also accept a previously computed dwellmap to speed up computation
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'ratemap',cell(1,1),'spatial_info_bsec',NaN,'spatial_info_bspike',NaN,'mutual_info',NaN,'sparsity',NaN,'spatial_coherence',NaN);
                sdatap = addToTable(sdatap,'place_fields',NaN,'field_area',cell(1,1),'field_centroids',cell(1,1),'field_weight_centroids',cell(1,1),'field_maj_lengths',cell(1,1),'field_min_lengths',cell(1,1),'field_orientations',cell(1,1),'field_snr',cell(1,1));
                sdatap = addToTable(sdatap,'grid_autocorrelation',cell(1,1),'grid_score',NaN,'grid_score2',NaN,'grid_wavelength',NaN,'grid_radius',NaN,'grid_orientation',NaN,'grid_mask',cell(1,1));
                sdatap = addToTable(sdatap,'over_dispersion_z',cell(1,1),'over_dispersion',NaN,'over_dispersion_r',NaN);
                          
                if (config.fild_asis || config.grid_asis || config.ovrd_asis) && numel(pspx)>minspikes % if we want to analyse place fields, grid cells or overdispersion
                    config.dwellmap = [];
                    if isfield(pdata.(part_now),'dwellmap') % if a previously computed dwellmap exists, use this to save time
                        config.dwellmap = pdata.(part_now).dwellmap;
                    end
                    % This function incorporates a number of different methods for generating firing rate maps and you can check the function for detailed descriptions of these
                    % my favourite is the method proposed by Leutgeb et al. (2005) Independent Codes for Spatial and Episodic Memory in Hippocampal Neuronal Ensembles
                    % https://doi.org/10.1126/science.1114037
                    % although the fastest method is 'nearest' as this is just a straight binning method followed by smoothing, so it is the default
                    [ratemap,dwellmap,~,config,mapdata] = mapDATA([ppox ppoy],[pspx pspy],config);

%% ########## %% ratemap analysis
                    % This function incorporates a number of different analyses and you can check the function for detailed descriptions of these
                    % I will point out that the spatial information content value provided by this function is calculated as in 
                    % Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
                    % https://doi.org/10.1002/(SICI)1098-1063(1996)6:2%3C149::AID-HIPO6%3E3.0.CO;2-K
                    % I note this because the method differs slightly between papers and I have found this one reflects the data the best
                    spatm = spatialMETRICS(ratemap,dwellmap);
                    
                    % accumulate data         
                    pdata.(part_now).dwellmap = single(dwellmap);
                    sdatap = addToTable(sdatap,'ratemap',{single(ratemap)},'spatial_info_bsec',spatm.spatial_information,'spatial_info_bspike',spatm.spatial_information_perspike,'mutual_info',spatm.mutual_info,'sparsity',spatm.sparsity,'spatial_coherence',spatm.spatial_coherence);

%% ########## %% place field analysis
                    % this analysis is pretty standard and is based on what I used before, a good citation would be:
                    % Park, Dvorak and Fenton (2011) Ensemble Place Codes in Hippocampus: CA1, CA3, and Dentate Gyrus Place Cells Have Multiple Place Fields in Large Environments
                    % "A place field was defined as any contiguous set of 9 or more pixels with greater that 0 AP/s firing rate that shared at least one side with another pixel in the field."
                    % https://dx.doi.org/10.1371%2Fjournal.pone.0022349
                    [fieldd,nfields] = getPFIELDS(ratemap,config.frcut,config.minfr,config.arcut);
                    % Although this function will take a minimum firing rate (i.e. replace the 0 AP/s) and I do tend to use 1Hz for that

                    % accumulate data
                    sdatap = addToTable(sdatap,'place_fields',nfields,'field_area',{fieldd.Area},'field_centroids',{fieldd.Centroid},'field_weight_centroids',{fieldd.WeightedCentroid},'field_maj_lengths',{fieldd.MajorAxisLength},'field_min_lengths',{fieldd.MinorAxisLength},'field_orientations',{fieldd.Orientation},'field_snr',{fieldd.signal_to_noise});
                    
%% ########## %% grid characteristic analysis
                    if config.grid_asis          
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
                        [~,gdata] = gridSCORE(automap);
                        
                        % accumulate data
                        sdatap = addToTable(sdatap,'grid_autocorrelation',{single(automap)},'grid_score',gdata.g_score,'grid_score2',gdata.g_score2,'grid_wavelength',gdata.wavelength,'grid_radius',gdata.radius,'grid_orientation',gdata.orientation,'grid_mask',{gdata.ring_mask});
                    end            

%% ########## %% Overdispersion analysis              
                    if config.ovrd_asis
                        % Taken from: Fenton et al. (2010) Attention-like modulation of hippocampus place cell discharge
                        % "In the first method, the entire session was divided into 5-sec intervals. For each interval we calculated the expected number of spikes, exp, as [equation in paper]
                        % where ri is the time-averaged rate at location i, and ti is the time spent in location i during the pass. Only intervals during which exp ? 5.0 AP were 
                        % used to calculate overdispersion since the overall firing rate of place cells is ~1.0 AP/sec.
                        % For each selected 5-sec interval, we then calculated z, the normalized standard deviation of obs, the observed number of spikes as [equation in paper]
                        % z measures the deviation of observed discharge from expected in standard deviation units. Overdispersion in turn is the variance of 
                        % the z distribution for a set of passes. The outcome of this calculation was found to be indistinguishable from the somewhat different 
                        % method previously used (Fenton and Muller, 1998)."
                        % https://doi.org/10.1523/JNEUROSCI.5576-09.2010
                        [overz,overd,overr] = getOVERDISPERSE([mapdata.poxnew mapdata.poynew],ppot,pspt,ratemap);
                        % essentially we compare the instantaneous firing of the cell at every position to the firing we would expect at that position according to the firing rate map
                        % the overdispersion z is the deviation of the expected from the observed, the overdispersion value is the standard deviation of this distribution
                
                        % accumulate data
                        sdatap = addToTable(sdatap,'over_dispersion_z',{overz},'over_dispersion',overd,'over_dispersion_r',overr);
                    end
                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% HD analyses    
                % HD firing rate maps are now computed by mapHD, like the above it will accept a precomputed
                % dwell time map to speed up computation
                % rayleigh vector values etc are computed on sum normalised HD firing rate maps to remove
                % any possible effects of firing rate (it's unclear if circ_stats care about the values)
                % If you recorded HD using two LEDs this should be used here and throughout the data
                % if not, HD is estimated using the animal's displacement
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'hd_spikemap',cell(1,1),'hd_ratemap',cell(1,1),'hd_rayleigh',NaN,'hd_max',NaN,'hd_mean',NaN,'hd_stdev',NaN);

                if config.head_asis && numel(pspx)>minspikes
                    % Load a head direction dwell map if one exists, or create it if not
                    if isfield(pdata.(part_now),'hd_dwellmap')
                        hd_dwell = pdata.(part_now).hd_dwellmap;
                        [~,~,hd_spikemap,hd_ratemap,r,mx,mn,sd] = mapHD(hd_dwell,ppoh,psph,config);
                    else
                        [~,hd_dwellmap,hd_spikemap,hd_ratemap,r,mx,mn,sd] = mapHD([],ppoh,psph,config);
                        pdata.(part_now).hd_dwellmap = single(hd_dwellmap);
                    end                

                    % accumulate data
                    sdatap = addToTable(sdatap,'hd_spikemap',{single(hd_spikemap(:))},'hd_ratemap',{single(hd_ratemap(:))},'hd_rayleigh',r,'hd_max',mx(1),'hd_mean',mn,'hd_stdev',sd);
                end
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Spike theta phase analyses   
                % By default any LFP analyses in klustest make use of the first LFP channel contained in
                % the mtint (i.e. LFP channel 1 in dacqUSB)
                % This can be changed at the top of klustest if you want, but it will crash if the LFP file does not
                % exist. Ideally the 'best' LFP channel would be used, where 'best' would be defined as the 
                % channel with the strongest theta or something. But for speed, we are just using the 1st one here
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'theta_phase_mean',NaN,'theta_phase_r',NaN,'theta_phase_max',NaN,'theta_phase_dist',cell(1,1));

                if config.spph_asis && numel(pspx)>minspikes
                    % bin the theta phase data
                    ai = reshape(deg2rad(-180:5:540),[],1); % doing this means we have bins symmetrical around zero
                    
                    % we must calculate the spike phase here precisely for every spike time, rather than indexing into the position thet phase
                    % this is for reasons of resolution (LFP and spikes are sampled at 48kHz positions are only sampled at 50Hz)
                    % This process is fairly well established:
                    % "To obtain a theta phase angle for each spike, LFPs were first bandpass filtered (fourth-order Chebyshev, r = 0.5, MATLAB 
                    % filter and filtfilt routines; 6–10 Hz) before a Hilbert transform was applied to obtain the instantaneous phase angle."
                    % from van der Meer and Redish (2011) Theta Phase Precession in Rat Ventral Striatum Links Place and Reward Information (https://doi.org/10.1523/JNEUROSCI.4869-10.2011)             
                    h = hilbert(lftheta);
                    phase = mod(angle(h),2*pi);  
                    pspp = interp1(lfpt(:),phase,'nearest');
                    spp2 = reshape([pspp; (pspp+2*pi)],[],1);
                    yi = histcounts(spp2,ai);
                    
                    % directional analyses on phase data
                    mx2p = ai(yi == max(yi)); % preferred angle (location of max frate)

                    % accumulate data
                    sdatap = addToTable(sdatap,'theta_phase_mean',circ_mean(pspp),'theta_phase_r',circ_r(pspp(:)),'theta_phase_max',mx2p(1),'theta_phase_dist',{yi});
                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Inter-spike interval analyses, theta and bursting analyses, refractory period analyses   
                % ISI analyses are quick and there shouldn't really be any reason to need to disable them,
                % however, it is assumed here that if you don't want refractory period or theta autocorrelogram analyses
                % ISI analyses are not required. The intrinsic theta fit can be seen in the klustest figure output
                % as a red line fitted to the theta spike autocorrelogram
                % This can also be used in combination with the frequency of global (LFP) theta as an indication of phase precession
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'isi_dist',cell(1,1),'isi_fdist',cell(1,1),'isi_fwhmx',NaN,'isi_half_width',NaN(1,4));
                sdatap = addToTable(sdatap,'burst_index',NaN,'burst_length_median',NaN,'burst_length_mean',NaN);
                sdatap = addToTable(sdatap,'intrinsic_theta_index',NaN,'intrinsic_theta_frequency',NaN,'intrinsic_theta_fit',NaN,'t500_spike_autocorr',cell(1,1),'t500_spike_autofit',cell(1,1));
                sdatap = addToTable(sdatap,'rpv_total',NaN,'rpv_proportion',NaN,'rpv_false_positive1',NaN,'rpv_false_positive2',NaN,'rpv_censored',NaN,'t25_spike_autocorr',cell(1,1));

                if config.sisi_asis && numel(pspx)>minspikes           
%% ########## %% inter-spike interval
                    % I don't remember exactly where I took this analysis from, but extracting the full width at half maximum is a fairly striahgtforward
                    % concept so I don't think it requires citation.
                    [~,idata,isis] = getISIhalfwidth(pspt);

                    % accumulate data
                    if isnan(idata.fwhmx) 
                        sdatap = addToTable(sdatap,'isi_dist',{NaN(1,101,'single')},'isi_fdist',{NaN(1,101,'single')},'isi_fwhmx',idata.fwhmx,'isi_half_width',NaN(1,4));
                    else
                        sdatap = addToTable(sdatap,'isi_dist',{single(interp1(idata.adist(:,1),idata.adist(:,2),0:0.5:50))},'isi_fdist',{single(interp1(idata.fdist(:,1),idata.fdist(:,2),0:0.5:50))},'isi_fwhmx',idata.fwhmx,'isi_half_width',[idata.half_max idata.fdist(idata.hwidth_ps,1)']);  
                    end

%% ########## %% burst index
                    % from Mizuseki et al. (2012) Activity Dynamics and Behavioral Correlates of CA3 and CA1 Hippocampal Pyramidal Neurons 
                    % "The spike-burst index was defined as the fraction of spikes with <6 ms ISIs (Harris et al., 2001)"
                    % https://doi.org/10.1002/hipo.22002
                    bindx = zeros(size(pspt));
                    bindx([isis; NaN] < 6 | [NaN; isis] < 6) = 1; % bindx is an index of all spikes sharing an isi less than 6ms
                    bindx = logical(bindx);
                    burst_index = sum(bindx) / numel(pspt);
                    sts = regionprops(bindx,'Area');

                    % accumulate data
                    sdatap = addToTable(sdatap,'burst_index',burst_index,'burst_length_median',nanmedian([sts.Area].'),'burst_length_mean',nanmean([sts.Area].'));

%% ########## %% Intrinsic theta    
                    % from van der Meer and Redish (2011) Theta Phase Precession in Rat Ventral Striatum Links Place and Reward Information
                    % "To quantify the degree and frequency of theta modulation in single cells, we used the method used by Royer et al. (2010). 
                    % First we computed the autocorrelogram of the cell, in 10 ms bins from -500 to +500 ms, normalized it to the maximum value
                    % between 100 and 150 ms (corresponding to theta modulation), and clipped all values above 1. Then we fit the following function [function in paper]
                    % where t is the autocorrelogram time lag, and a-c, w, and t1–2 were fit using the fminsearch optimization function in MATLAB. A measure
                    % of theta modulation strength, the “theta index,” was defined as a/b, which intuitively corresponds to the ratio of the sine fit relative to
                    % the baseline in the autocorrelogram."
                    % https://doi.org/10.1523/JNEUROSCI.4869-10.2011
                    [thi,thf,thr,~,c500,f500] = getTHETAintrinsic(pspt); % we also extract the frequency of theta and the goodness of the theta fit

                    % accumulate
                    sdatap = addToTable(sdatap,'intrinsic_theta_index',thi,'intrinsic_theta_frequency',thf,'intrinsic_theta_fit',thr,'t500_spike_autocorr',{single(c500(:))'},'t500_spike_autofit',{single(f500(:))'});
                    
%% ########## %% Refractory period analyses           
                    % refractory period violation analyses
                    % from Navratilova et al. (2016) Grids from bands, or bands from grids? An examination of the effects of single unit contamination on grid cell firing fields.
                    % although they took these analyses from pre-existing papers
                    % "There are few reliable methods for estimating the contamination of a unit isolated from tetrode recordings (see DISCUSSION). The only method that does not rely 
                    % on the same measures that are used for spike sorting is to check the number of spikes that occur within the refractory period of another spike..."
                    % https://doi.org/10.1152/jn.00699.2015
                    [nrpv,prpv,fp1,fp2,censored,~,c25] = getRPVcount(pspt,isis,part_duration); % we extract the number of spikes in the refractory period, what proportion of the total this is etc                    

                    % accumulate data
                    sdatap = addToTable(sdatap,'rpv_total',nrpv,'rpv_proportion',prpv,'rpv_false_positive1',fp1,'rpv_false_positive2',fp2,'rpv_censored',censored,'t25_spike_autocorr',{single(c25(:))'});
                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Speed modulation analysis   
                % This analysis is taken from Kropff et al (mEC speed cell paper)
                % It is unlikely to be of great use to HPC people, although many place cells are speed modulated
                % Related to this, you can also compute the relationship between theta power and running speed
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'speed_score',NaN,'speed_slope',NaN,'speed_y_intercept',NaN,'speed_frate_curve',cell(1,1));
                
                if config.sped_asis && numel(pspx)>minspikes
                    [svals,stime,sscore,sslope,sintcpt,crve,~] = getSPEEDmod(ppot,ppov,pspt);

                    % accumulate data
                    pdata.(part_now).speed_dwell_time = single(interp1(svals,stime,0:1:50,'linear'));     
                    sdatap = addToTable(sdatap,'speed_score',sscore,'speed_slope',sslope,'speed_y_intercept',sintcpt,'speed_frate_curve',{single(interp1(svals,crve,0:1:50,'linear'))});
                end                
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Cell type/identity  
                % I use stuff like this to help categorise cells, so I have included it here in case others find it useful
                % Realistically, for final analyses you will probably want to test SI and grid score etc against a shuffle
                % rather than a hard cutoff. But for exploration (which klustest is geared towards) it is fine for a guide
                % It will work best when the cutoffs in getCELLTYPE2 are best suited to your needs

                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'cell_type',cell(1,1),'cell_type_int',NaN,'cell_type_bin',cell(1,1));
                
                if config.cell_asis && config.wave_asis && config.fild_asis && config.grid_asis && config.head_asis && numel(pspx)>minspikes % we can only do cell typing with all this info available
                    [btype,ctype,ttype] = getCELLTYPE2(sdatap.frate,sdatap.waveform_params(1,2),sdatap.spatial_info_bsec,sdatap.grid_score,sdatap.hd_rayleigh);                

                    % accumulate data
                    sdatap = addToTable(sdatap,'cell_type',{ttype},'cell_type_int',ctype,'cell_type_bin',{btype});
                end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Part figure and finish loop
                % figPART is the main figure, it creates a summary for each cluster in each part, using the info in mtint, pdata and sdatap
                % I have tried to optimise it as far as possible
                
                if numel(pspx)>minspikes
                    if run_figPART
                        figPART(mtint,pdata,sdatap,fig_vis); % overall figure with ratemaps etc
                    end
                end

                % concatenate sdatap (this part data) into sdatac (this cluster's data)
                if ~isempty(sdatap)
                    sdatac = [sdatac;sdatap];
                end

%% ########## %% Add data to bdata table (behaviour data)
                if isempty(bdata) || ~any(bdata.partn==pp) % only add behaviour data if this part doesn't exist in bdata yet
                    % accumulate, these variables don't need to be preallocated as they should always be filled
                    bdatap.rat = {pdata.rat};
                    bdatap.date = str2double(pdata.date);
                    bdatap.partn = pp;
                    bdatap.part = {part_now};                
                    bdatap.duration = part_duration;
                    bdatap.directory = {pwd};
                    bdatap = addToTable(bdatap,'positions',{single([ppox ppoy])},'speed',{single(ppov)},'HD',{single(ppoh)},'dwellmap',{single(dwellmap)},'HDdwellmap',{single(pdata.(part_now).hd_dwellmap)});

                    % concatenate bdatap (this part data) into bdata (this cluster's data)
                    bdata = [bdata;bdatap];
                end
            end % this ends the parts loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Cluster figures and finish loop            
            if run_figCLUS
                if nparts > 1 % if there is more than one part (otherwise this figure is not useful)
                    figCLUS(mtint,pdata,sdatac,fig_vis); % overall figure with ratemaps etc
                end 
            end

            % concatenate sdatac (this cluster's data) into sdata (data for all clusters)
            if ~isempty(sdatac)
                sdata = [sdata;sdatac];
            end
            loopout = looper(loopout);

        end % this ends the cluster loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Electrode figures and finish loop  

        %% Cluster cross-correlation figure
        if run_figCROSS
            figfile = [pwd '\klustest\' pdata.combined_name '\figures\'];
            [~,~,~] = mkdir(figfile);                
            figCROSS(sdata,'tet',tet,'figfile',figfile,'fig_vis','off')
        end
    % 
    %     %% Cluster space figure
    %     if run_figCSPACE
    %         [~,~,~] = mkdir([pwd '\klustest\' sdata.combined_name '\figures\']);
    %         figfile = [pwd '\klustest\' sdata.combined_name '\figures\' sdata.combined_name '_E' num2str(tet) '_cluster_space.png'];
    %         figCSPACE(sdata.(tetstr).fetdata,figfile,fig_vis); % overall figure with ratemaps etc
    %     end
    end % this ends the electrode loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the data and finish up
    save([pwd '\klustest\' config.cname '\' config.cname '_sdata.mat'],'sdata'); % save session data
    save([pwd '\klustest\' config.cname '\' config.cname '_pdata.mat'],'pdata'); % save session data
    save([pwd '\klustest\' config.cname '\' config.cname '_bdata.mat'],'bdata'); % save session data

    % finish up
    toc1 = toc/60;
    disp(sprintf('klustest has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',[pwd '\klustest\' config.cname '\figures'],' &'');">','figures folder','</a>'])
    disp('-------------------------------------------------------------------------------------');

end % ends the main klustest function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A stupid subfunction for preallocating a table
function tin = addToTable(tin,varargin)
    for i=1:2:length(varargin)
        tin.(varargin{i}) = varargin{i+1}; % table in (variable name) = variable value to assign
    end
end % ends addToTable function







