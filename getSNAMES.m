function [snames,cname,nsess,fnames] = getSNAMES(get_all)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function just gathers some session names by asking the user to select .set files
%   it outputs some useful forms of this info
%   [snames,cname,nsess] = getSNAMES
%
%%%%%%%% Outputs
%   snames = cell array of selected file names, filename with no extension or path
%   cname = a new filename generated using all the input snames concatenated
%   nsess = the number of sessions/set files selected
%   fnames = a cell array of full file names
%
%%%%%%%% Comments 
%   24/11/16 created so this functionality can be shared between 2D and 3D functions
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial settings
if ~exist('get_all','var') || isempty(get_all)
    get_all = 0;
end % if ~exist('get_all','var') || ismepty(get_all)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find sessions
if get_all
    snames = dir('*.set');
    if ~numel(snames)
        error('ERROR: no detected .set files, check working directory matches data directory... exiting');
    end
    snames = {snames(:).name};
    snames = snames(:);
else
    snames = uipickfiles('FilterSpec','*.set','Output','cell','Prompt','Select the sessions to analyse...');
end % if get_all
nsess = numel(snames);

% sort out parameters for later and reduce .set file names
fnames = {};
if length(snames) == 1
    nnow = snames{1};
    [~,nme,~] = fileparts(nnow);
    cname = nme;  
    snames{1} = nme;    
    fnames{1} = [pwd '\' nme];
else
    for ff = 1:length(snames)
        nnow = snames{ff};
        [~,nme,~] = fileparts(nnow);
        snames{ff} = nme;
        fnames{ff} = [pwd '\' nme];
        if ff == 1 % the first filename
            cname = nme;
        elseif ff == length(snames) % the last filename
            cname = [cname '_' nme];
        else % middle filenames
            cname = [cname '_' nme '_'];
        end % if ff == 1 % the first filename
    end % for ff = 1:length(snames)
end % if length(fnames) == 1






