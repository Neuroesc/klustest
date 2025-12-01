function [fwhmx,idata,isi] = getISIhalfwidth(spt,twin,isi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes a spike train and calculates full width at half maximum of spike ISI
%   It also returns full width at tenth maximum and an estimation of the standard deviation of the FWHM
%   [fwhmx,stdv,fwtmx] = forNING(spt,twin)
%
%%%%%%%% Inputs
%   spt = spike times in seconds, if omitted function will generate a fake spike dataset based on sum of gaussians
%   twin = (default = 60) time window over which to compute half width, around 50ms seems to work well
%
%%%%%%%% Outputs
%   fwhmx = full width at half maximum of spike ISI
%   stdv = 1 standarddeviation of spike ISI, estimated from fwhmx
%   fwtmx = full width at one tenth maximum of spike ISI
%
%%%%%%%% Comments
%   31/10/17 created as a toy example for Ningyu
%   14/11/17 advanced to a subfunction for Klustest
%   17/11/17 added catch for when there is only 1 spike (no interspike intervals)
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
% if ~exist('spt','var') || isempty(spt) % if no spikes are provided, make a fake spike dataset
%     nspikes = 500;
%     ds = [normrnd(10,10,nspikes,1) normrnd(112,20,nspikes,1) normrnd(225,30,nspikes,1) normrnd(350,40,nspikes,1)];
%     spt = cumsum(ds); % turn into spike train times
%     spt = spt/1000; % convert to seconds
% end

if ~exist('twin','var') || isempty(twin) || all(isnan(twin(:)))
    twin = 60; % time window for calculation and figures
end
bin_size = 1;

spt = spt(:).*1000; % convert to ms

idata = struct;
idata.no_window = 0;
idata.adist = [reshape(movmean(0:bin_size:twin,2,'EndPoints','discard'),[],1) zeros(numel(movmean(0:bin_size:twin,2,'EndPoints','discard')),1)];
idata.fdist = [reshape(0:bin_size:twin,[],1) zeros(numel(0:bin_size:twin),1)];
idata.half_max = NaN;
idata.hwidth_ps = [NaN NaN NaN];
idata.fwhmx = NaN;
fwhmx = NaN;
idata.stdev = NaN;
idata.burst_index_6ms = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate FWHM
% calculate inter spike intervals
if ~exist('isi','var') || isempty(isi) || all(isnan(isi(:)))
    isi = diff(spt);
end
if isempty(isi)
    isi = NaN;
    return
end
idata.burst_index_6ms = sum(isi<6) ./ numel(isi);
idata.burst_index_25ms = sum(isi<25) ./ numel(isi);

% ISI histogram
edg = 0:bin_size:twin;
f = histcounts(isi,edg);
cents = movmean(edg,2,'EndPoints','discard');
idata.adist = [cents(:) f(:)];

% fit distribution to data
pd = fitdist(isi,'Kernel','Kernel','epanechnikov','Bandwidth',4); % get a smoothed density estimate of the data (otherwise its too irregular)
x2 = 0:bin_size:twin; % new time intervals over which to calculate distribution
ipdf = pdf(pd,x2); % get new distribution
ipdf = ipdf ./ nansum(ipdf) .* 0.5; % normalise to unit integral
idata.fdist = [x2(:) ipdf(:)];

if ~any(ipdf) % if there are no isis within the pdf window it will be all zeros
    return
end

% FWHM full width at half maximum
[totalMax,indx0] = findpeaks(ipdf,'NPeaks',1); % find first peak (maximum)
halfMax = totalMax / 2; % find the half max value
hindx1 = find(ipdf(indx0:end) <= halfMax,1,'first')+(indx0-1); % find where the data first drops below half the max, after the max

if isempty(hindx1) || isempty(indx0)
    return
end

hindx2 = find(ipdf(1:indx0) <= halfMax,1,'last'); % find where the data is last below half the max, before the max
if isempty(hindx2)
    hindx2 = 1;
end
fwhmx = x2(hindx1) - x2(hindx2); % the full width at half maximum is the distance between these two points
stdv = fwhmx./(2*sqrt(2*log(2))); % standard deviation from hlaf width

% accumulate data
idata.half_max = halfMax;
idata.hwidth_ps = [indx0 hindx2 hindx1];
idata.fwhmx = fwhmx;
idata.stdev = stdv;

%% Generate figures if necessary
if 0 % set to 0 to switch figures off
    % plot stuff
    figure('visible','on','position',[100 100 1200 800]);
    subplot(2,2,1)
    bar(cents,f,1,'k')
    hold on
    bar(-cents,f,1,'k')
    xlabel('Time (ms)')
    ylabel('Count')
    title('500ms ISI histogram')    
    
    subplot(2,2,2)
    x3 = 0:0.5:500; % new time intervals over which to calculate distribution
    ipdf3 = pdf(pd,x3); % get new distribution
    ipdf3 = ipdf3 * 0.5; % multiply by bin width so it sums to 1    
    area(x3,ipdf3,0,'FaceColor',[0.5 0.5 0.5])
    hold on
    area(-x3,ipdf3,0,'FaceColor',[0.5 0.5 0.5])
    xlabel('Time (ms)')
    ylabel('Probability')    
    title('500ms KDE ISI');
    
    subplot(2,2,3)
    edg = 0:0.5:twin;
    f = histcounts(isi,edg);
    cents = histcents(edg);
    bar(cents,f,1,'k')
    hold on
    bar(-cents,f,1,'k')
    xlabel('Time (ms)')
    ylabel('Count')
    title(sprintf('%.fms ISI histogram',twin))

    subplot(2,2,4)
    area(x2,ipdf,0,'FaceColor',[0.5 0.5 0.5])
    hold on
    area(-x2,ipdf,0,'FaceColor',[0.5 0.5 0.5])
    xlabel('Time (ms)')
    ylabel('Probability')
    title(sprintf('%.fms KDE ISI',twin))

    hold on
    ax = gca;
    line([x2(hindx2) x2(hindx1)],[halfMax halfMax],'Color','k','LineWidth',1.5)
    line([x2(indx0) x2(indx0)],ax.YLim,'Color','k','LineWidth',1)
    text(x2(hindx1),halfMax,sprintf('FWHM: %.1f ms',fwhmx),'FontSize',8)
    text(x2(indx0),halfMax*2,sprintf('Max: %.3f',halfMax*2),'FontSize',8) 
end



















