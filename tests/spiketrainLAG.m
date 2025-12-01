function [spikesOUT,lagOUT] = spiketrainLAG(spikesIN,tmax,lagIN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes a spike train (given as time points in s) and lags it by a certain amount (+ve or -ve)
%   or it will lag it by a random amount
%   [spikesOUT,lagOUT] = spiketrainLAG(spikesIN,tmax,lagIN)
%
%%%%%%%% Inputs
%   spikesIN = the original spike times
%   tmax = the actual maximum time of the recording (defaults to the max time in tin)
%   lagIN = the amount of lag wanted in seconds (defaults to a random value +/- tmax-1
%
%%%%%%%% Outputs
%   spikesOUT = the output lagged spike times
%   lagOUT = the lag which was used
%
%%%%%%%%%%% log
%   23/04/16 created for calculating shuffles spatial information content
%   24/06/16 overhauled for Ele
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('spikesIN','var') || isempty(spikesIN)
    error('No spikes presented... exiting')
end % if ~exist('in','var') || ismepty(in)
if ~exist('tmax','var') || isempty(tmax)
    tmax = nanmax(tin);
end % if ~exist('in','var') || ismepty(in)
if ~exist('lin','var') || isempty(lagIN)
    lag = randi(tmax-1,[1 1]);
else
    lag = lagIN;
end % if ~exist('in','var') || ismepty(in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate lag
% Get the number of spikes
spiken = length(spikesIN);

% Triplicate differences between spike times
trip_spikes = diff([0 spikesIN]); % running difference index
trip_spikes = repmat(trip_spikes,1,3); % triplicate difference index

% Offset spike train by desired amount
trip_spikes_offset = trip_spikes + lagIN;

% Extract offset spike train
trip_mid = trip_spikes_offset(spiken+1:spiken*2); % cut out middle section of triplicate difference index after offset
spikesOUT = cumsum(trip_mid); % take cumulative sum to regain spike times that should now be offset
lagOUT = lag;















