function [part_indx,part_names,part_times] = partSESSION(mthd,mtint,part_names,part_ignore)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes some input information and returns parameters specifying how to partition
%   a recording session
%   [part_indx,part_names,part_times] = partSESSION(mthd)
%
%%%%%%%% Inputs
%   mthd = the partitioning style, 1 to partition based on inp file, 2 to partition into seperate sessions, 3 to partition everything into one file (whole session)
%   mtint = the mtint produced by klustest etc
%
%%%%%%%% Outputs
%   part_indx = an index of which parts to include (a partition may be formed from other partitions, i.e. if we want a 'whole session' output as well as parts)
%   part_names = the names of the partitions, may be given as an input or maybe generated here automatically
%   part_times = the start and end time of each partition
%
%%%%%%%% Comments
%   24/11/16 created so this functionality could be shared between 2D and 3D functions
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('part_names','var') || isempty(part_names)
    part_names = [];
end % if ~exist('in','var') || ismepty(in)
if ~exist('part_ignore','var') || isempty(part_ignore)
    part_ignore = ones(length(part_names));
end % if ~exist('in','var') || ismepty(in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine partitions
% sort out array of start and end times
if mthd == 1 % if we want to divide the output based on digital inputs
    % sort out vector of goal identifiers
    goal_array = mtint.inps.goal; % load goal values
    part_indx{1} = find(goal_array == 1); % find all these values in goal array
    part_indx{2} = find(goal_array == 2); % find all these values in goal array							                                        
    part_indx{3} = find(goal_array == 3); % find all these values in goal array								                                        
    part_indx{4} = find(goal_array == 4); % find all these values in goal array								                                        
    part_indx{5} = cell2mat(part_indx); % find all these values in goal array
    part_names = {'part1' 'part2' 'part3' 'part4' 'Entire'}; 

    part_times = mtint.inps.tstamps; % timestamp values from file
    part_times = part_times(:); % make sure it is a column
    if mod(numel(part_times),2)
        error('ERROR: there is an odd number of digital inputs... exiting')
    end % if mod(x,2)
    part_times = [part_times(1:2:end,:) part_times(2:2:end)]; % place start and end times side by side    

elseif mthd == 2 % if we want to divide the output based on input sessions
    dvals = mtint.pos.trial_duration;
    dvals = [0; dvals(:)];
    dvals = cumsum(dvals);
    dvals = [dvals circshift(dvals,-1)];
    dvals(end,2) = sum(dvals(:,1));
    dvals = dvals(1:end-1,:);

    part_times = dvals;
    part_indx = num2cell(1:length(dvals(:,1)));
    if isempty(part_names)
        for ss = 1:numel(snames)
            part_names{ss} = ['s_' snames{ss}];
        end % for ss = 1:numel(snames)
    end % if isempty(part_names)
    part_times = part_times(part_ignore == 1,:);
    part_names = part_names(part_ignore == 1);
    part_indx = num2cell(1:numel(part_names));
    
elseif mthd == 3
    part_indx{1} = 1;
    part_times = [0 sum(mtint.pos.trial_duration)];
    part_names{1} = ['s_' cname]; % append session name with an s_ to avoid structure fields beginning with a number (i.e. an error)
    
else
    error('ERROR: unknown method input to partSESSIONS... exiting')

end % if part_sessions






























