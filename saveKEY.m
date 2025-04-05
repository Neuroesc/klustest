function [count,timestamps,type,value,text_value,intervals] = saveKEY(iname,meth,config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function saves or loads a .key file (a text file version of an .inp file)
%   [count,timestamps,type,value] = saveKEY(iname,meth)
%
%%%%%%%% Inputs
%   iname = the name of the .inp file (or key file, it doesn't really matter)
%   meth = (default = 1) 1 to save a .key file for a corresponding .inp file, 2 to load the contents of an existing .key file
%
%%%%%%%% Outputs
%   count = the number of inputs
%   timestamps = the time values for these
%   type = the type of each
%   value = the value of each (uint32 integer)
%   text_value = the value of each (actual string or keyboard value)
%
%%%%%%%% Comments
%   19/05/16 created 
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('iname','var') || isempty(iname) || all(isnan(iname))
    error('ERROR: a filename must be provided to saveKEY... exiting')
end

if ~exist('meth','var') || isempty(meth) || all(isnan(meth))
    meth = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process .key file
% read .inp file
[~,nme,~] = fileparts(iname);
inpname = [pwd '\' nme '.inp'];
keyname = [pwd '\' nme '.txt'];
    
if meth == 1 % if we want to save a .key file
    if ~exist(inpname,'file')
        error(sprintf('ERROR: .inp file %s not found... unable to load',inpname))
    end
    [count,timestamps,type,value] = read_key(inpname);

    % process the value data to fix errors
    value = cellfun(@transpose,value,'UniformOutput',0); 
    value = cellfun(@num2str,value,'UniformOutput',0);
    value = cellfun(@strtrim,value,'UniformOutput',0);
    value = regexprep(value,{'\W',' '},'');
    value(cellfun(@isempty,value)) = {666}; % some keys are not recognised by uint32, replace these with NaNs, usually this happens when we hit a weird key by accident like ` or ¬
    value = value';
    %value = cellfun(@uint32,value,'UniformOutput',0);
    
    % find interval keys so we can highlight them
    intervals = cell(length(timestamps(:)),1);
    if exist('config','var')
        idx = ismember(value,config.interval_keys);
        intervals(idx) = {'<--'};
    end

    % write a .key file
    x = table;
    x.count = (1:length(timestamps(:)))';
    x.timestamps = cellstr(num2str(timestamps'));
    x.type = uint32(cell2mat(type(:)));
    x.value = cell2mat(cellfun(@uint32,value(:),'UniformOutput',0));
    x.text_value = cell2mat(value(:));
    x.intervals = intervals;
    %x.text_value = sprintf('%c',x.value(:));
    writetable(x,keyname,'FileType','text','Delimiter','tab');

    % read .key file back (to ensure function outputs are always the same)
    xout = readtable(keyname); 
    count = xout.count(:);
    timestamps = xout.timestamps(:);
    type = xout.type(:);
    value = xout.value(:);
    text_value = xout.text_value(:);
    if ~iscell(text_value)
        text_value2 = cell(length(text_value),1);
        for i = 1:length(text_value)
            text_value2{i,1} = text_value(i,:);
        end
        text_value = text_value2;
    end
    text_value = cellfun(@num2str,text_value,'UniformOutput',0);
    intervals = xout.intervals(:);

elseif meth == 2 % if we want to load from a .key file
    % read .key file
    xout = readtable(keyname,'FileType','text','Delimiter','tab'); 
    count = xout.count(:);
    timestamps = xout.timestamps(:);
    type = xout.type(:);
    value = xout.value(:);
    text_value = xout.text_value(:);
    if ~iscell(text_value)
        text_value2 = cell(length(text_value),1);
        for i = 1:length(text_value)
            text_value2{i,1} = text_value(i,:);
        end
        text_value = text_value2;
    end
    text_value = cellfun(@num2str,text_value,'UniformOutput',0);
    intervals = xout.intervals(:);
    
end




















































