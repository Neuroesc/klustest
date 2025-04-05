function [tets,snames,cutname] = getTRODES(config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getTRODES  find the names of the sessions contributing to this cutfile and the active tetrodes
% Given the combined name of a session this function checks which tetrode has a cut file, opens it and extracts the name of the 
% files (sessions) that contributed towards it. In this process we also get a list of which tetrodes have a .cut file (i.e. were
% cluster cut).
%
% USAGE:
%         [tets,mtets,snames,cutname] = getTRODES(cname,tin)
%
% INPUT:
%         cname - combined name of the session, without an extension, default is 'kwiktint' because thats the default output for kwiktint
%         tin - (optional) tetrodes to check for, default is 1:16
%
% OUTPUT:
%    tets - list of tetrodes with cut files
%    snames - cell array of session names
%    cutname - the name of the cutfile that was used
%
% EXAMPLES:
%
% See also: klustest getcut setdiff

% HISTORY:
% version 1.0.0, Release 24/08/16 Initial release
% version 1.0.1, Release 04/04/18 Commenting and formatting for klustest update
% version 2.0.0, Release 04/04/18 Changed to also search for active tetrodes (combined getTRODES and getCNAMES)
% version 2.0.1, Release 04/04/18 Removed list of inactive tetrodes, this can be generated using setdiff, reduces vargin to just 1
% version 2.0.2, Release 19/04/19 Updated to remove possible duplicates of tetrodes
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
% find which tetrodes are present
cnames = dir([config.cname '*.cut*']); % list all cut files corresponding to the given combined name
clist = {cnames.name}.'; % get their names
cnums = cellfun(@(x) regexp(x,'[_]\d+[.]','match'),clist); % extract just the number part of the name
tets = cellfun(@(x) str2double(regexp(x,'\d+','match')),cnums); % convert to an integer, these numbers are the existing tetrodes (cut files)
tets = unique(tets); % sometimes there are backed up .cut files, so tetrodes are duplicated, if these are not removed klustest will run multiple times on each tetrode

% open the first existing tetrode and extract the session info
cutname = [config.cname '_' num2str(tets(1)) '.cut']; % the cut file name
[~,etext] = getcut(cutname); % get data from cutfile, specifically the line which identifies the parent sessions
idx = strfind(etext,': '); % should contain two values
flist = textscan(etext(idx(1)+2:idx(2)-8),'%s','delimiter',','); % get the file name parts of the string
snames = flist{1,1}; % extract the relevant cell
    
    
    
    
    
    
    
    
    
    
    
    







