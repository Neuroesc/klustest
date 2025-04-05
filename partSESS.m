function part_config = partSESS(part_config,mname,posn,skip_trial_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes some input information and returns parameters specifying how to partition a recording session
%   part_config = partSESS(part_config,mtint)
%
%%%%%%%% Inputs
%   part_config = part_config structure with fields specifying the name, method, intervals and times etc for each part
%
%%%%%%%% Outputs
%   part_config = part_config structure with added part times and inputs (if needed)
%
%%%%%%%% Comments
%   24/11/16 created so this functionality could be shared between 2D and 3D functions
%   01/04/17 changed to a cell array format
%   06/04/17 changed to a structure format
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if isfield(part_config,'FDATA')
    nparts = length(fieldnames(part_config))-1;
else
    error('ERROR: FDATA field missing from part_config, this should contain recording times... exiting');
end % if isfield(part_config,'FDATA')

if ~exist('skip_trial_fig','var')
    skip_trial_fig = 0;
end

%% Get the recording starts and ends
dvals = part_config.FDATA.recording_times; % get the duration of each recording session
dvals = [0; dvals(:)]; % add a zero value
dvals = cumsum(dvals); % calculate cumulative total, this gives the start/end of each session within the combined files
dvals = [dvals circshift(dvals,-1)]; % we want to go from a single vector of end times to a matrix where each row contains a start and end time
dvals = dvals(1:end-1,:); % we don't want the last row because this just represents the end of the last trial

%% Fix keypresses for the whole recording, if they exist
fnames = fieldnames(part_config);
mths = NaN(nparts,1);
for pp = 1:nparts
    part_method = part_config.(fnames{pp}).method;
    mths(pp) = part_method;
end
if any(mths==3)
    [stimes,etimes] = figKEYS(mname,part_config,skip_trial_fig,posn);
    atimes = [stimes(:) etimes(:)]; % all time start-end pairs in rows
end

%% Get the times for each part
frem = []; % index of fields to remove
fnames = fieldnames(part_config);
for pp = 1:nparts
    part_method = part_config.(fnames{pp}).method;
    
    if part_method == 1 % if the method is whole session
        part_config.(fnames{pp}).times = [min(dvals(:)), max(dvals(:))]; % start time and end time equals start and end of whole file

    elseif part_method == 2 % if the method is recordings
        recvec = part_config.(fnames{pp}).intervals; % vector specifying which recordings should be included
        if max(recvec) > numel(dvals(:,1))
            disp(sprintf('\t\tERROR: requested recording (%d) is greater than the number of recordings (%d)... skipping',max(recvec),numel(dvals(:,1))));
            frem = [frem pp];
            continue
        end 
        part_config.(fnames{pp}).times = dvals(recvec(:),:); % start and end times are the start and end times of the requested recording sessions

    elseif part_method == 3 % if the method is digital inputs
        % determine which trials to include
        inpvec = part_config.(fnames{pp}).intervals; % vector specifying which input intervals should be included
        if isinf(inpvec) % if inpvec is infinite just use all trials
            inplog = ones(size(stimes));
        elseif max(inpvec) > numel(stimes)
            error('ERROR: requested input interval (%d) is greater than the number of input pairs (%d)... exiting',max(inpvec),numel(stimes));
        else % normal handling of trials
            inplog = zeros(size(stimes));
            inplog(inpvec) = 1;
        end 
        
        % add data to part_config 
        part_config.(fnames{pp}).times = atimes(logical(inplog),:); % start and end times are the start and end times of the requested recording sessions 
    end 
end 
for rr = 1:length(frem)
    part_config = rmfield(part_config,fnames{frem(rr)});
end





































