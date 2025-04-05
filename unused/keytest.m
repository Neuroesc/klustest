% function keytest(in,in2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION  short descr.
% long descr.
%
% USAGE:
%         [out] = template(in,in2)
%
% INPUT:
%         in - input 1
%         in2 - input 2
%
% OUTPUT:
%    p - output
%
% EXAMPLES:
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with input variables
% deal with the other variables
inps = {'tetrodes','cname'};
vals = {'1:16','''kwiktint'''};
for ff = 1:length(inps)
    if ~exist(inps{ff},'var')
        eval([inps{ff} '=' vals{ff} ';']);
    end
end

% get the function name
stk = dbstack;
function_name = stk.name;

% initialise settings
maintain_mtint = 0; % set to 1 to save/load mtint in the base workspace, this saves time when running the function mutliple times (for instance in debugging) but should otherwise be set to 0
mtint_override = 1; % set this to 1 to force the production of a new mtint file (do this if you reclustered in tint for instance)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get session data
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


% check all inp files are present
inpPresent = zeros(1,length(snames));
for i = 1:length(snames)
    inpPresent(i) = exist([pwd '\' snames{i} '.inp'],'file');
end % for i = 1:length(snames)

% read .inp files and accumulate their data
inps = struct;
inps.count = 0; 
inps.tstamps = []; 
inps.value = [];
inps.index = [];
inps.spec = [];
mtint.header = getDACQDATAHEADERS(snames,'.set');

ttime = 0;
filepath = [pwd '\'];
[fileStruct,tetsAvailable] = listDACQFILES(filepath,snames);
if all(inpPresent)
    for pp = 1:numel(fileStruct)      
%         if exist([snames{pp} '.key'],'file') % if a .key file already exists
%             [count,timestamps,spec,value] = saveKEY([snames{pp} '.inp'],2); % load it
%         else
            [count,timestamps,spec,value] = saveKEY([snames{pp} '.inp'],1); % if not then make one     
%         end           

        inps.count = inps.count + count;
        timestamps = timestamps + ttime;
        inps.tstamps = [inps.tstamps; single(timestamps(:))];
        inps.value = [inps.value; value(:)];
        inps.index = [inps.index; ones(size(timestamps(:))).*pp];
        inps.spec = [inps.spec; spec(:)];
        timenow = mtint.header(pp).duration;
        ttime = ttime + timenow;
    end 
end










keyboard








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% retrieve keys
disp(sprintf('\t...retrieving data from mtint'));
load(mname,'inps');
load(mname,'pos');
vals = inps.value;
tims = inps.tstamps;
pox = pos.xy_pixels(:,1); % position x
poy = pos.xy_pixels(:,2); % position y
pot = pos.ts; % position times

% run through each cell and do bullshit conversion that we need
disp(sprintf('\t...converting to strings'));
v2 = cell(size(vals));
for ii = 1:length(vals)
    vn = vals{ii};
    if ~isstring(vn)
        vn = num2str(vn);
    end
    vn = sprintf('%s',vn);
    vn = strrep(vn,' ','');
    v2{ii} = vn;
end










