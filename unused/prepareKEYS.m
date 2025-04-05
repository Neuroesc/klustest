function prepareKEYS(varargin)
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
%           'tetrodes'          -   (default = 1:16) Vector, the tetrodes to run on (i.e. [1 2 3 4 5 6] would run on tetrodes 1 to 6)
%
%           'clusters'          -   (default = 0) Vector, the clusters to run on, set this to 0 to run on all clusters
%
%           'rname'             -   (default = first 3 characters of parent directory) String, the rat name/number; this will be used in the sdata structure so its important to give this
%
%           'cname'             -   (default = 'kwiktint') String, function will automatically look for klustakwiked files named this
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
% version 01.0.0, Release 15/01/20

%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2016 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INPUT ARGUMENTS CHECK
%% Prepare default settings 
    def_cname               = 'kwiktint';   
    % I can extract the rat name from the file path because I know my folder structure, if this is not the case for you, change def_rname to what you prefer
    % the next two lines make rname equal the first 3 characters of the directory containing the current directory
    pname = pwd; sindx=strfind(pname,'\'); 
    def_rname = pname(sindx(end-1)+1:sindx(end-1)+3);    
    
%% Parse inputs
    p = inputParser;
    addOptional(p,'cname',def_cname,@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    addOptional(p,'rname',def_rname,@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    parse(p,varargin{:});

%% Retrieve parameters 
    config = p.Results;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INITIAL SETTINGS / INPUTS
    config.interval_keys    = {'s','e'}; % what keypresses were used to delineate trial starts and ends {start,end} these can be integers or characters i.e. {1,2} or {'s','e'}, {'1','2'} is the same as {1,2}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% KEY FILES
    disp(sprintf('Assessing data...'))
    [~,snames,~] = getTRODES(config.cname); 




        manageKEYS(snames,2,{'s','e'});





















