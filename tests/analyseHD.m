function [hdfr,ahb,shb,spike_hd,rvec,pfd_mean,pfd_std,pfd_max,max_fr] = analyseHD(hd,bins,srate,stimes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A short function to analyse head direction data
% hd = head direction data for the whole session, in (+ve) degrees
% bins = the angular bins to use for hd (default = 0:30:360)
% srate = sampling rate of the position data in ms
% stime = spike times of the neuron to be analysed
%
% 01/03/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(bins,1) < size(bins,2)
    bins = bins';
end % if size(bins,1) < size(bins,2)
if size(hd,1) < size(hd,2)
    hd = hd';
end % if size(bins,1) < size(bins,2)

post = (0:length(hd)) * (srate/1000);
post = post';
post = round(post*100);
stimes = round(stimes*100);
spike_hd = hd(ismember(post,stimes));

ahb = histc(hd,bins) * (srate/1000); % to get binned hd sampling
shb = histc(spike_hd,bins); % to get binned spikes in hd

hdfr = shb ./ ahb; % to get binned hd firing rate

rvec = circ_r(bins,hdfr);
pfd_mean = circ_mean(bins,hdfr);
[i,j] = max(hdfr);
max_fr = i;
pfd_max = bins(j);
pfd_std = circ_std(bins,hdfr);

