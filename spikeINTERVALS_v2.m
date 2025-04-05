function [tisi,sisi,hisi] = spikeINTERVALS_v2(spt,wsiz,bsiz,spt2,smoo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function calculates a spike autocorrelation, given spike times and a window length
%   If two spike trains are given, it can also calculate the cross-correlogram.
%   In the cross-correlogram, a peak before zero would mean spikes in the second spike train occur 
%   frequently right after the spikes in the first train
%   [tisi,sisi,hisi] = spikeINTERVALS_v2(spt,wsiz,bsiz,spt2)
%
%%%%%%%% Inputs
%   spt = spike times (s)
%   wsiz = window size (ms) default is 400ms
%   bsiz = binsize (ms) default is 1ms
%   spt2 = optional second spike train, then a cross-correlogram is produced intead of an auto-correlogram
%
%%%%%%%% Outputs
%   tisi = time (ms) of correlogram bins
%   sisi = spike probability
%   hisi = spike histogram
%
%%%%%%%% Comments
%   22/08/16 created from spk_acorr2 because I also want to be able to cross-correlate different cells
%   22/08/16 added cross-correlation method, vectorised some operations
%   06/04/16 further vectorised, replaced clumsy windows with bsxfun
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
% if no secondary spikes were given
if ~exist('spt2','var') || isempty(spt2) || all(isnan(spt2)) || ~any(spt2) % if no spikes were given
    spt2 = spt;
end

% if the binsize was not specified
if ~exist('bsiz','var') || isempty(bsiz) 
    bsiz = 0.5;
end

% if the window width was not specified
if ~exist('wsiz','var') || isempty(wsiz) 
    wsiz = 400;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike histogram
spt = spt(:);
spt2 = spt2(:);
nspikes = size(spt2,1); % the number of spikes in the second train

hedge = -wsiz:bsiz:wsiz; % histogram bin edges
tisi = mean([hedge(1:end-1); hedge(2:end)]); % histogram bin centres
hisi = zeros(size(tisi));
sisi = zeros(size(tisi));
if ~exist('spt','var') || isempty(spt) || all(isnan(spt)) || ~any(spt) % if no spikes were given
    return
end

if max([size(spt,1) size(spt2,1)]) > 2^12 % if there are so many spikes we cannot make an square array
    for ii = 1:numel(spt2) % for every spike in train 2
        sdiffs = spt-spt2(ii);
        sdiffs(~sdiffs) = NaN;
        sdiffs = sdiffs * 1000;
        
        hisit = histcounts(sdiffs(:),hedge);
        hisi = hisi + hisit; % accumulate spike counts
    end
else
    sdiff = bsxfun(@minus,repmat(spt,1,nspikes),spt2');
    sdiff(~sdiff) = NaN;
    sdiff = sdiff * 1000; % convert to ms

    hedge = -wsiz:bsiz:wsiz;
    hisi = histcounts(sdiff(:),hedge);
end
sisi = hisi ./ nspikes;















