function [snames,cutname1] = getCNAMES(cname,tin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getCNAMES  find the names of the sessions contributing to this cutfile
% Given the name of a cut file this function checks which tetrode has a cut file, opens it and extracts the name of the 
% files (sessions) that contributed towards it
%
% USAGE:
%         [snames,cutname1] = getCNAMES(cname,tin)
%
% INPUT:
%         cname - name of the .cut file, without the extension, default is 'kwiktint' because thats the default output for kwiktint
%         tin - (optional) tetrodes to check for, default is 1:16
%
% OUTPUT:
%    snames - cell array of session names
%    cutname1 - the name of the cutfile used
%
% EXAMPLES:
%
% See also: klustest setdiff fileparts

% HISTORY:
% version 1.0.0, Release 24/08/16 Initial release
% version 1.0.1, Release 04/04/18 Commenting and formatting for klustest update
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with input variables
inps = {'cname','tin'};
vals = {'kwiktint','1:16'};
reqd = [0 0];
for ff = 1:length(inps)
    if ~exist(inps{ff},'var')
        if reqd(ff)
            error('ERROR: vargin %s missing... exiting',inps{ff});            
        end        
        eval([inps{ff} '=' vals{ff} ';']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
cnames = dir([cname '*.cut*']); % list all cut files corresponding to the given combined name
clist = {cnames.name}.'; % get their names
cnums = cellfun(@(x) regexp(x,'[_]\d+[.]','match'),clist); % extract just the number part of the name
cnums = cellfun(@(x) str2double(regexp(x,'\d+','match')),cnums); % extract just the number part of the name









for tt = 1:length(tin) % for every hypothetical tetrode
    cutname1 = [cname '_' num2str(tin(tt)) '.cut']; % the cut file name
    if any(exist(cutname1,'file')) % see if it exists
        break % if it exists use this filename
    end
end
[~,etext] = getcut(cutname1); % get data from cutfile, specifically the line which identifies the parent sessions
idx = strfind(etext,': '); % should contain two values
flist = textscan(etext(idx(1)+2:idx(2)-8),'%s','delimiter',','); % get the file name parts of the string
snames = flist{1,1}; % extract the relevant cell
disp(sprintf('\t...working on sessions: %s',strjoin(snames,', ')));
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    







