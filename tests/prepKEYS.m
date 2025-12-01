function prepKEYS(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
%prepKEYS function to load digital keypresses and visualise them 
%    This function is designed to run before klustest etc have had their chance to analyse the data
%    It is for quick assessment and correction of digital keypresses
%
% USAGE:
%           prepKEYS() process files with default settings
%                 prepKEYS will look for .cut files named [cname] (default is kwiktint), it will open the first one and extract 
%                 details about which sessions were used to make that .cut file. It will then open all of the required pos files
%                 and key files if they exist or .inp files if they do not. It will then process the digitial keypresse and plot
%                 some summary plots that can be used to verify/check keypresses. There is an option to load the key files, edit
%                 them, save them and reload the data with these new values.
%
%           prepKEYS(Name,Value,...) process with Name-Value pairs used to control aspects of the function output
%
%           Parameters include:
%
%           'mtint'              -   Structure, mtint structure produced by klustest
%
%           'cname'              -   ['kwiktint'] String, function will automatically look for klustakwiked files named this
%
%           'intervals'          -   [{'s','e'}] what keypresses were used to delineate trial starts and ends {start,end} these can be integers or characters i.e. {1,2} or {'s','e'}, {'1','2'} is the same as {1,2}
%
%           'ignore_cut'         -   [false] Logical scalar, true means the function will ask you to specify which .inp files to use, false means the function will open a .cut file
%                                            with a filename that matches 'cname' and extract the session names from that, then open the .inp files matching those session names
%
%           'overwrite'          -   [2] Scalar, 1 will cause new key files to be created, overwriting old ones, 2 will use old key files if they exist
%
% EXAMPLES:
%
%           % run function using default values
%           prepKEYS()
%
%           % run function using default values, but overwrite existing key files
%           prepKEYS('overwrite',1)
%
%           % run function using default values, but specify a different filename to use for .cut files
%           klustest('cname','data_file_name')
%
% See also: klustest manageKEYS

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
    def_intervals           = {'s','e'};   
    def_ignore_cut          = false;   
    def_mtint               = struct();   
    def_overwrite           = 2;   

%% Parse inputs
    p = inputParser;
    addParameter(p,'cname',def_cname,@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    addParameter(p,'intervals',def_intervals,@(x) ~isempty(x) && iscell(x) && length(x)==2); 
    addParameter(p,'ignore_cut',def_ignore_cut,@(x) ~isempty(x) && numel(x)==1);     
    addParameter(p,'mtint',def_mtint,@(x) ~isempty(x) && isstruct(x));   
    addParameter(p,'overwrite',def_overwrite,@(x) ~isempty(x) && numel(x)==1);         
    parse(p,varargin{:});

%% Retrieve parameters 
    config = p.Results;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% KEY FILES
    disp(sprintf('Identifying sessions...'))
    if ~isempty(fieldnames(config.mtint))
        manageKEYS(config.mtint,config.overwrite,config.intervals,0);        

    else
        flist = dir('*.cut*');
        cnames = {flist.name}';
        if isempty(cnames) || config.ignore_cut % if there are no cut files or if we want to choose our own sessions
            fnames = uigetfile('*.inp*','Select .inp files','MultiSelect','on');
            snames = fnames';
            for pp = 1:length(snames)
                [~,n,~] = fileparts(snames{pp});
                snames{pp} = n;
            end
            disp(sprintf('\t...working on sessions: %s',strjoin(snames,', ')));    
        else
            tetrodes = 1:32;
            for tt = 1:length(tetrodes) % for every hypothetical tetrode
                cutname1 = [config.cname '_' num2str(tetrodes(tt)) '.cut']; % the cut file name
                if any(exist(cutname1,'file')) % see if it exists
                    break % if it exists use this filename
                end
            end
            disp(sprintf('\t...reading %s',cutname1));
            [~,etext] = getcut(cutname1); % get data from cutfile, specifically the line which identifies the parent sessions
            idx = strfind(etext,': '); % should contain two values
            flist = textscan(etext(idx(1)+2:idx(2)-8),'%s','delimiter',','); % get the file name parts of the string
            snames = flist{1,1}; % extract the relevant cell
            disp(sprintf('\t...working on sessions: %s',strjoin(snames,', ')));    
        end
        manageKEYS(snames,config.overwrite,config.intervals,0);
        
    end




















