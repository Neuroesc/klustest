function [ts,ch1,ch2,ch3,ch4] = get_spikes(filename,ffast)
%
%   [ts,ch1,ch2,ch3,ch4] = getspikes(filename)
%   
%   Copyright (C) 2004 Sturla Molden 
%   Centre for the Biology of Memory
%   NTNU
%
[spikes,spikeparam] = importspikes(filename);
ts = cell2mat({spikes(:).timestamp1})';
ch1 = NaN;
ch2 = NaN;
ch3 = NaN;
ch4 = NaN;
if ~exist('ffast','var') || isempty(ffast) || ~ffast
    nspk = spikeparam.num_spikes;
    spikelen = spikeparam.samples_per_spike;
    ch1 = reshape(cell2mat({spikes(:).waveform1}),spikelen,nspk)';
    ch2 = reshape(cell2mat({spikes(:).waveform2}),spikelen,nspk)';
    ch3 = reshape(cell2mat({spikes(:).waveform3}),spikelen,nspk)';
    ch4 = reshape(cell2mat({spikes(:).waveform4}),spikelen,nspk)';
end % if isempty(ffast) || ~exist(ffast,'var') || ~ffast
















