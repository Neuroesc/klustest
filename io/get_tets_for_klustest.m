function [tets,snames,data_dirs] = get_tets_for_klustest(dataformat,config)
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
    switch dataformat
        case {'kwikcut','neuralynx'}
            % find possible directories
            possible_dnames = dir(pwd);
            possible_dnames = {possible_dnames([possible_dnames.isdir]).name};
            possible_dnames = possible_dnames(~ismember(possible_dnames,{'.','..'}));
            
            % find which directories have spike data
            data_dirs = [];
            snames = [];
            for dd = 1:length(possible_dnames)
                f1 = dir([pwd '\' possible_dnames{dd} '\*.ntt']); % neuralynx tetrode file
                f2 = dir([pwd '\' possible_dnames{dd} '\*.nst']); % neuralynx stereotrode file
            
                if ~isempty(f1) || ~isempty(f2)
                    data_dirs = [data_dirs; {[pwd '\' possible_dnames{dd}]}];
                    snames =  [snames; possible_dnames(dd)];
                end
            end
        
            % find tetrodes in each sub directory
            p_tet = 1:1024;
            for ff = 1:length(data_dirs) % for every data directory
                d = dir([data_dirs{ff} '\*.ntt']);
                d([d.bytes]==16384) = []; % remove files that only contain a header
                fnames_ntt = {d.name};
                tets_found = cell2mat(cellfun(@str2num,replace(fnames_ntt,{'TT','.ntt','.NTT',},''),'UniformOutput',false));
                p_tet = intersect(p_tet,tets_found); % if an electrode is missing in one directory it will be discarded here
            end
    
            % find stereotrodes in each sub directory
            p_ste = 1:1024;
            for ff = 1:length(data_dirs) % for every data directory
                d = dir([data_dirs{ff} '\*.nst']);
                d([d.bytes]==16384) = []; % remove files that only contain a header
                fnames_ntt = {d.name};
                tets_found = cell2mat(cellfun(@str2num,replace(fnames_ntt,{'ST','.nst','.NST',},''),'UniformOutput',false));
                p_ste = intersect(p_ste,tets_found); % if an electrode is missing in one directory it will be discarded here
            end
            tets = [ [p_tet(:) zeros(size(p_tet(:)))]; [p_ste(:) ones(size(p_ste(:)))] ];

        case {'kwiktint'}
            % find which tetrodes are present
            cnames = dir(['*.cut*']); % list all cut files corresponding to the given combined name
            clist = {cnames.name}.'; % get their names
            cnums = cellfun(@(x) regexp(x,'[_]\d+[.]','match'),clist); % extract just the number part of the name
            tets = cellfun(@(x) str2double(regexp(x,'\d+','match')),cnums); % convert to an integer, these numbers are the existing tetrodes (cut files)
            tets = unique(tets); % sometimes there are backed up .cut files, so tetrodes are duplicated, if these are not removed klustest will run multiple times on each tetrode
            
            % find which tetrodes we want to (and can) analyse            
            tets = intersect(config.tetrodes,tets);
            tets(:,2) = zeros(size(tets(:,1)));

            % open the first existing tetrode and extract the session info
            cutname = [config.cname '_' num2str(tets(1)) '.cut']; % the cut file name
            [~,etext] = getcut(cutname); % get data from cutfile, specifically the line which identifies the parent sessions
            idx = strfind(etext,': '); % should contain two values
            flist = textscan(etext(idx(1)+2:idx(2)-8),'%s','delimiter',','); % get the file name parts of the string
            snames = flist{1,1}; % extract the relevant cell

            % find where the original data files are stored 
            data_dirs = strcat(repmat(pwd,size(snames,1),1),repmat('\',size(snames,1),1),snames); 

        case {'phy'}
            % find where the original files are stored
            cnames = dir(pwd); % list all cut files corresponding to the given combined name
            data_dirs = {cnames.name}.'; % get their names    
            data_dirs = data_dirs(~ismember(data_dirs,{'klustest','.','..'}) & [cnames.isdir]');

            % find which tetrodes are present            
            cnames = dir([data_dirs{1} '\' data_dirs{1} '.mountainsort\output*']); % list all cut files corresponding to the given combined name
            clist = {cnames.name}.'; % get their names
            cnums = cellfun(@(x) regexp(x,'\d+','match'),clist,'UniformOutput',false); % extract just the number part of the name
            tets = cellfun(@(x) str2double(x),cnums); % convert to an integer, these numbers are the existing tetrodes (cut files)
            tets = unique(tets); % sometimes there are backed up .cut files, so tetrodes are duplicated, if these are not removed klustest will run multiple times on each tetrode
            
            % find which tetrodes we want to (and can) analyse            
            tets = intersect(config.tetrodes,tets);
            cutname = 'phy';

            % session info (need to add)
            snames = data_dirs; % extract the relevant cell
    end
    
    
    
    
    
    
    
    
    
    
    







