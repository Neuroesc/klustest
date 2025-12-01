function part_config = getPARTconfig(fname,part_config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function generates a figure showing the cluster space of a tetrode
%   klusspace(tet,mtint,sdata,fig_vis)
%
%%%%%%%% Inputs
%   fname = the filename, either where we want to save part_config or where we want to load it from
%   method = 1 for saving part_config to fname, 2 to load part_config from fname
%   part_config = the part_config cell array
%
%%%%%%%% Comments
%   04/04/17 created to contain this loading crap
%   06/04/17 updated to newer structure format, removed 3D exception as its not needed now
%   06/04/17 fixed issue, where function would not read last line of saved part_config
%   10/05/17 added more detailed handling of method 3
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('fname','var') || isempty(fname)
    error('ERROR: getPARTconfig requires at least 1 input... exiting');
end % if ~exist('fname','var') || isempty(fname) || ~exist('method','var') || isempty(method) 

if nargin == 1 % if only a filename was given, we probably want to load the thing
    method = 2;
else
    method = 1;
end % if nargin == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run through file
if method == 1 % save the structure as a .txt file
    fname2 = [fname '_' datestr(now,30)];
    if exist(fname,'file')
        disp(sprintf('\tcopying %s to %s...',fname,fname2));
        copyfile(fname,fname2,'f');
    end % if exist(fname,'file')

    fid = fopen(fname,'w'); % open it and assign it an ID
    part_names = fieldnames(part_config); % get the upper field names of the structure - these should be the part names
    for pp = 1:length(part_names) % for every part name
        pn = part_names{pp}; % get the current part name
        field_names = fieldnames(part_config.(pn)); % get the field names under this part name

        fprintf(fid,'session_heading=%s\r',pn); % print the part name after 'name='
        for ff = 1:length(field_names) % for every one of the sub fields
            fn = field_names{ff}; % get its name
            fv = num2str(part_config.(pn).(field_names{ff})); % get its value
            fprintf(fid,'%s=%s\r',fn,fv); % print the value after the name followed by an '='
        end % for ff = 1:length(field_names)
    end % for pp = 1:length(part_names)
    fclose(fid); % close the text file

elseif method == 2 % load the structure from a .txt file
    fid = fopen(fname,'r'); % open it and assign it an ID
    datout = textscan(fid,'%s','Delimiter','\r'); % scan in the text
    datout = datout{1,1}; % get the text part we want
    fclose(fid); % close the text file
    nindx = strfind(datout,'session_heading'); % find lines containing the 'name' data
    nindx(cellfun('isempty',nindx)) = {NaN}; % fill empty index values with NaN
    nindx = find(cell2mat(nindx)==1); % find where the non-NaNs are
    nindx = [nindx; length(datout)+1]; % add the last line as an index, we need this later
    part_names = {}; % empty an array for data
    part_config = struct; % empty an array for data
    for nn = 1:length(nindx)-1 % for every line containing 'session_heading'
        dnow = datout{nindx(nn)}; % get it
        dnow = strsplit(dnow,'='); % split based on the '='
        part_names{nn} = dnow{1,2}; % the part name is the second piece after the '='
        for ll = nindx(nn)+1:nindx(nn+1)-1 % for every line between this one and the next one (or the end of the file)
            lnow = datout{ll}; % get it
            lnow = strsplit(lnow,'='); % split based on the '='
            if size(lnow,2)==1 % if there was nothing after the '='
                lnow{1,2} = []; % set it to empty
            else
                charval = lnow{1,2}(isstrprop(lnow{1,2},'alpha')); % if there was something after the '=', find how many letter it contains
                if ~length(charval) && ~strcmp(lnow{1,1},'out_name') % if there are no alphanumeric characters and this is not the outname
                    lnow{1,2} = str2num(lnow{1,2}); % convert the contents to digits
                end % if ~charval
            end % if size(lnow,2)==1
            part_config.(part_names{nn}).(lnow{1,1}) = lnow{1,2}; % place the data in a structure in the correct field}
        end % for ll = nindx(nn):nindx(nn+1)-1
    end % for nn = 1:length(nindx)
end % if method == 1







