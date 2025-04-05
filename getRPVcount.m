function [nrpv,prpv,fp1,fp2,censored,t25,c25] = getRPVcount(pspt,isis,dur,t25,c25)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getRPVcount  calculate refractory period violations
% count the number and find the proportion of spikes which fall within a cells refractory period
% can be used as an estimate of cluster quality
%
% USAGE:
%         [nrpv,prpv,fp1,fp2,censored,t25,c25] = getRPVcount(pspt,isis,dur,t25,c25)
%
% INPUT:
%         pspt - vector of spike times
%         isis - interspike interval, can be calculated using getISIhalfwidth
%         dur - duration of the session
%         t25 - time lag values for a spike autocorrelogram
%         c25 - probability values for a spike autocorrelogram
%
% OUTPUT:
%    nrpv - number of refractory period violations
%    prpv - proportion  of total spikes which are refractory period violations
%    fp1 - false positive rate 1
%    fp2 - false positive rate 2
%    censored - the number of spikes missed because of the system lockout
%    t25 - time lag values for the spike autocorrelogram used
%    c25 - probability values for the spike autocorrelogram used
%
% See also: getISIhalfwidth klustest

% HISTORY:
% version 1.0.0, Release 09/04/19 Initial release, contained here to clean up klustest
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2019 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with input variables
narginchk(3,5);

% refractory period settings
tau_r = 2; % length of refractory period in ms
tau_c = 0.75; % lockout period of recording system in ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
    if ~exist('t25','var') || isempty(t25) || ~exist('c25','var') || isempty(c25) % if a spike autocorrelogram was not provided, make one
        [c25,t25] = spike_auto(pspt,'win_size',25,'bin_size',0.5);
    end

    % refractory period violation analyses
    tmt = tau_r - tau_c; % max time in which we can observe a noise spike
    nspikes = numel(pspt); % number of cluster spikes
    frate = nspikes ./ dur; % firing rate of cluster in Hz
    nrpv = sum(isis <= tau_r); % number of refractory period violations
    prpv = nrpv ./ nspikes; % refractory period violations as a proportion of all spikes

    % Method 1 taken from Dan Hill's UltraMegaSort
    % Spikes occurring under the minimum ISI are by definition false positives, their total rate relative to the rate of the neuron overall gives the estimated false positive rate.
    vtime = 2 * nspikes * tmt; % total time available for violations - there is an available window of length (tau_r - tau_c) after each spike.
    vrate = nrpv/vtime;
    fp1 = vrate/frate;
    if fp1 > 1
        fp1 = 1; % happens sometimes and is actually an error in the assumptions of the analysis
    end

    % Method 2 using equation from Hill et al. (2011) Quality metrics to accompany spike sorting of extracellular signals
    % This method might return imaginary results for large numbers of violations
    k = 2 * tmt * (nspikes^2);
    rts = roots([k -k nrpv*dur]);
    fp2 = min(rts);

    censored = (tau_c/1000) * nspikes / dur; % estimate the fraction of spikes not detected because of the system lockout






























