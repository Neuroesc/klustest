function [tisi,sisi,hisi] = spikeINTERVALS(spt1,wind,bsiz,spt2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function calculates a spike autocorrelation, given spike times and a window length
%   If two spike trains are given, it can also calculate the cross-correlogram. This is calculated as the
%   probability of spikes in train2 occuring at spike times in train1
%   [tisi,sisi] = spikeINTERVALS(spt1,bsiz,wind,spt2)
%
%%%%%%%% Inputs
%   spt1 = spike train 1, given in ms time stamps
%   wind = the window length in ms, this is the width of the resulting cross- or auto-correlogram
%   spt2 = (optional) a second spike train, if we want to calculate the cross-coorelogram between this and the first one
%
%%%%%%%% Outputs
%   tisi = time (ms) of correlogram bins
%   sisi = spike probability
%   hisi = spike histogram
%
%%%%%%%% Comments
%   22/08/16 created from spk_acorr2 because I also want to be able to cross-correlate different cells
%   22/08/16 added cross-correlation method, vectorised some operations
%   © Roddy Grieves: rmgrieves@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('spt1','var') || isempty(spt1) % if no spikes were given
    error('ERROR: no spike times given to spikeINTERVALS... exiting')
end % if ~exist('in','var') || ismepty(in)

% if the binsize was not specified
if ~exist('bsiz','var') || isempty(bsiz) 
    bsiz = 1;
end % if ~exist('bsiz','var') || isempty(bsiz) % if the binsize was not specified
bsiz = bsiz*1000; % convert to ms

% if the window width was not specified
if ~exist('wind','var') || isempty(wind) 
    wind = 500;
end % if ~exist('wind','var') || isempty(wind) % if the window width was not specified

% if we have been given two spike trains
if exist('spt2','var') && ~isempty(spt2)
    mode = 2;
else
    mode = 1;
end % if length(spt(1,:)) == 2 % if we have been given two spike trains

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sort out the data first
n_spikes = numel(spt1); % total number of spikes
max_time = max(spt1); % the last spike time
n_bins = ceil(max_time * bsiz); % the number of bins to use when binning spike times

spt_round = ceil(spt1 * bsiz); % rounded spike times in ms
spk_train = zeros(n_bins,1); % an empty array to contain a logical index of spikes
spk_sorted = tabulate(spt_round); % bin the spike times into ms
spk_train(spk_sorted(:,1)) = spk_sorted(:,2); % place these values in the spike train vector

if mode == 2 % if we are cross-correlating
    spt_round2 = ceil(spt2 * bsiz); % rounded spike times in ms
    spk_train2 = zeros(n_bins,1); % an empty array to contain a logical index of spikes
    spk_sorted2 = tabulate(spt_round2); % bin the spike times into ms
    spk_train2(spk_sorted2(:,1)) = spk_sorted2(:,2); % place these values in the spike train vector    
else
    spk_train2 = spk_train;
end % if mode == 2 % if we are cross-correlating

%% For every spike, overlay a window of zeros and sum how many spikes fall into each of the windows bins
counts = zeros(wind+1,1);
counted = 0;
for b = 1:length(spt_round) % for every spike
    bin = spt_round(b); % get the rounded spike time which is also the bin index

    if ((bin > wind/2+1) && (bin < length(spk_train)-wind/2+1)) % as long as the window does not overlap the beginning or end of the session
        counts = counts + spk_train2(bin-wind/2:bin+wind/2); % sum the spikes observed in this window with the counts vector
        counted = counted + sum(spk_train2(bin-wind/2:bin+wind/2)) - spk_train2(bin); % sum the total number of spikes added to the counts vector so far, minus the spikes the window was centred on
    end % if ((bin > wind/2+1) && (bin < length(train)-wind/2+1))

end % for b = 1:length(spt_round) % for every spike

counts(wind/2+1) = 0; % remove the central count, which is the spikes the window was centred on
if max(counts) == 0 && counted == 0
    counted = 1;
end % if max(counts) == 0 && counted == 0
hisi = counts;
sisi = counts ./ counted;
tisi = -wind/2:1:wind/2;
    







