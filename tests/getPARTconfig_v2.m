function part_config = getPARTconfig_v2(fname,part_config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function completes a part_config structure with time data and also saves it as a text file
%   In another method it can also load the part_config from a txt file (handy for editing)
%   part_config = getPARTconfig_v2(fname,part_config)
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
%   19/05/17 changed a lot to take advantage of eval, this means the text file contains an exact readout of the structure, which makes editing it easier
%   19/05/17 changed to version 2
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


    fid = fopen(fname,'wt'); % open it and assign it an ID
    part_names = fieldnames(part_config); % get the upper field names of the structure - these should be the part names
    for pp = 1:length(part_names) % for every part name
        % the current part name
        pn = part_names{pp};

        % get all sub-field names
        field_names = fieldnames(part_config.(pn)); % also get the field names under this part name

        for ff = 1:length(field_names) % for every one of the sub fields
            % the current sub-field name
            fn = field_names{ff};

            % make strings for saving the data, the fieldname and it's value
            field_name_full = ['part_config.' pn '.' fn];
            
            if iscell(part_config.(pn).(fn))
                field_value = '{';
                for cc = 1:length(part_config.(pn).(fn))
                    if cc == 1
                        field_value = [field_value mat2str(part_config.(pn).(fn){cc})];                        
                    else
                        field_value = [field_value ',' mat2str(part_config.(pn).(fn){cc})];
                    end
                end
                field_value = [field_value '}'];
            else
                field_value = mat2str(part_config.(pn).(fn));
            end
            save_string = [field_name_full ' = ' field_value];

            % write to text file
            fprintf(fid,'%s\n',save_string); % print the value after the name followed by an '='
        end % for ff = 1:length(field_names)
    end % for pp = 1:length(part_names)
    fclose(fid); % close the text file

elseif method == 2 % load the structure from a .txt file
    fid = fopen(fname); % open it and assign it an ID
    part_config = struct; % make sure an empty part_config is ready

    % run through text file and recreate part_config
    tline = fgetl(fid); % get the first line of text file
    while ischar(tline) % while there is still text
        eval([tline ';']) % directly evaluate the line of text as a command, it should be something like 'part_config.part.method = 4'
        tline = fgetl(fid); % get the next line
    end
    fclose(fid);

end % if method == 1











































