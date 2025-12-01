function [jtimes,tshuff,sos] = shuffleSPIKES(duration,itimes,mint,pos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function randomly shuffles spike trains by a random time value (+ve or -ve)
%   If the position data (x,y,t) for the session are also given 
%   jtimes = shuffleSPIKES(duration,itimes,mint)
%
%%%%%%%% Inputs
%   duration = (s) the maximum length of the recording session
%   itimes = (s) the spike times
%   mint = (default = 20s) the minimum time we want to shuffle the train by
%   pos (optional) = the position data [x y t] coordinates
%
%%%%%%%% Outputs
%   jtimes = the new, shuffled spike train
%   tshuff = the time the train was shuffled by
%   sos (if pos is given) = the new spike [x y t] coordinates that correspond to the observed data, calculated using a nearest neighbour approach
%
%%%%%%%% Comments
%   10/08/16 created 
%   11/08/16 actually got it working
%   06/11/19 use shiftSTRAIN instead
%   © Roddy Grieves: rmgrieves@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
% if the important inputs are missing, display an error
if ~exist('duration','var') || isempty(duration) || ~exist('itimes','var') || isempty(itimes)
    error('ERROR: missing inputs to shuffleSPIKES... exiting');
end % if ~exist('duration','var') || isempty(duration) || ~exist('itimes','var') || isempty(itimes)

% if a minimum time for shuffling was not given, use 20s
if ~exist('mint','var') || isempty(mint)
    mint = 20;
end % if ~exist('mint','var') || ismepty(mint)
maxt = duration-0.02; % the maximum time is also the duration of the recording (-0.02 for safety)

% check the orientation of the itimes vector
if length(itimes(:,1)) > length(itimes(1,:))
    itimes = itimes';
end % if length(itimes(:,1)) < length(itimes(1,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform shuffle
% generate random time interval
rshuff = mint + (maxt-mint) .* rand(1,1); % generate random time shifts between mint and maxt (these are all positive)
rnegs = round(rand(1,1));
rnegs(~rnegs) = -1; % we need a vector of random 1s and -1s because we want to shift the trains forward and backwards
rshuff = rshuff .* rnegs;
tshuff = round(rshuff);

% move the spike train
jtimes = itimes + rshuff;

% check for negative spikes
zindx = find(jtimes <= 0);
jtimes(zindx) = duration - abs(jtimes(zindx));

% check for spikes over the duration
pindx = find(jtimes >= duration);
jtimes(pindx) = jtimes(pindx) - duration;

% re-sort the spike times
jtimes = sort(jtimes,'ascend');
jtimes = jtimes';

%% Work out the new spike positions
if exist('pos','var')
    pox = double(pos(:,1));
    poy = double(pos(:,2));
    pot = double(pos(:,3));
    n_n_indx = knnsearch(pot,double(jtimes)); % find the position times that correspond to the new spike times (nearest neighbour)
    spt_new = pot(n_n_indx);
    spx_new = pox(n_n_indx); % also extract the corresponding position data
    spy_new = poy(n_n_indx);
    sos = [spx_new spy_new spt_new];
    
else
    sos = NaN(1,3);
end % if exist('pos','var')


