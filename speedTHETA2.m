% function speedTHETA(cname)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes a klustest sdata structure as input and 
%   [out] = template(in,in2)
%
%%%%%%%% Inputs
%   in = 
%
%%%%%%%% Outputs
%   out =
%
%%%%%%%% Comments
%   00/00/00 created 
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
% if ~exist('in','var') || isempty(in) || all(isnan(in))
%     in = 
% end
close all

eeg_chann = 1; % which eeg channel to use
ffreqs = [6 12];

% spectrogram settings
pad = 0; % do we want to zero pad the FFT
err = [1 0.95]; % error values for computing FFT
tapers = [3 5]; % spectrogram tapers
wind_v = 1; % spectrogram window size (s)
step_v = 0.25; % spectrogram step duration (s)
range_v = [0.1 300]; % frequency range of spectrogram

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Theta speed calculations
% cname is the name to use when seeking cluster cut data (i.e. the name of your cut file minus the extension).
% specified in the next if section is the value of cname to use
% I have set up kwiktint so it outputs stuff named 'kwiktint' by default in most cases, so this should be the default here. 
% If you changed this in kwiktint you should change it here though.
% If you want cname to have a different value for a one off, pass it in as an input argument or initialise cname here instead:
% cname = 'b1'; % comment/edit/uncomment me when needed
% either way these next lines will be skipped and your input value will be used
if ~exist('cname','var') || isempty(cname) || all(isnan(cname))
    cname = 'kwiktint'; % (default = 'kwiktint') please don't edit me!
end
load(['klustest3\' cname '\' cname '_sdata.mat']); % load sdata
load(['klustest3\' cname '\' cname '_mtint.mat']); % load mtint (for eeg)

% part_config stuff
part_config = sdata.part_config;
nparts = length(fieldnames(part_config))-1;
part_names = fieldnames(part_config);

%% For each part, calculate the theta speed ratio
for pp = 1:nparts % for every partition
    part_now = part_names{pp}; % the name of the current part
    part_now3 = [part_now '_3d'];
    part_times = part_config.(part_now).times; % the time pairs (intervals) corresponding to this part

    % get the running speed for this part
    pov = sdata.(part_now3).pov; 
    pot = sdata.(part_now3).pot;

    % get the LFP data for the whole session (time should match pot though)
    lfp = double(mtint.lfp(eeg_chann).lfp(:,1));
    Fs = mtint.lfp(eeg_chann).Fs(1,1);
    lt = (0:length(lfp)-1)'/Fs; % make a vector for time


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % get the LFP for just this part
    part_times = part_config.(part_now).times; % the time pairs (intervals) corresponding to this part
    plfp = [];
    for ii = 1:length(part_times(:,1)) % for every pair of time values (interval) associated with this part  
        i_time = part_times(ii,:); % get the start and end time
        pindx = find(lt > i_time(1) & lt < i_time(2)); % find the position data falling into this interval
        plfp = [plfp; lfp(pindx,1)];
    end    

    % find PSD
    [pwindow,freqs] = pwelch(plfp(:),Fs,Fs/2,[],Fs);
    sdata.(part_now3).lfp.psd = [freqs(:) pwindow(:)];
    sdata.(part_now3).lfp.psd_theta_power = nanmean(pwindow(freqs >= 6 & freqs <= 12));
    [~,mp] = nanmax(pwindow(freqs >= 6 & freqs <= 12));
    freq2 = freqs(freqs >= 6 & freqs <= 12);
    freq2 = freq2(mp);
    sdata.(part_now3).lfp.psd_theta_frequency = freq2; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate figure
fig_spec = figure('visible','on','Position',[100, 100, 1024, 800]);
fig_hor = 3; % how many plots wide should it be
fig_ver = 1+nparts; % how many plots tall should it be
fspac = 0.01; % the spacing around the plots, on all sides
fpadd = 0.03; % the spacing around the plots, on all sides, this takes more space than fspac though
fmarg = 0.03; % the margins around the plots, at the edge of the figure
set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white

subaxis(fig_ver,fig_hor,1:3,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg); 
cla
imagesc(sp1a(1:30,:));      
axis xy
xlabel('Time (s)') % label x axis
ylabel('Frequency (Hz)') % label y axis
c = caxis;
caxis([0 c(2)]);

for pp = 1:nparts % for every partition
    part_now = part_names{pp}; % the name of the current part
    part_now3 = [part_now '_3d'];
    
    subaxis(fig_ver,fig_hor,fig_hor*pp+1,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);    
    freqs = sdata.(part_now3).lfp.psd(:,1);
    psdv = sdata.(part_now3).lfp.psd(:,2);
    plot(freqs,psdv);
    ax = gca;
    ax.XLim = [0 50];
    tf = sdata.(part_now3).lfp.psd_theta_frequency;
    tp = sdata.(part_now3).lfp.psd_theta_power;    
    title(sprintf('theta power = %.2f, frequency = %.2f',tp,tf),'FontSize',8)    
    ylabel('Power (dB)');
    xlabel('Frequency (Hz)');
    hold on
    line([6;6],ax.YLim,'Color','k');
    line([12;12],ax.YLim,'Color','k');
    line([tf;tf],ax.YLim,'Color','g');
    
    subaxis(fig_ver,fig_hor,fig_hor*pp+2,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);  

    
    subaxis(fig_ver,fig_hor,fig_hor*pp+3,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);    

end

% save figure
id = 'speedTHETA_analysis';
print(fig_spec,'-dpng','-r150',[pwd '\klustest3\' sdata.combined_name '\' id '.png'])
% close(fig_spec);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finish up
% Save the session data structure
info = whos('sdata');
siz = info.bytes / 1000000;
disp(sprintf('\t...saving %.1fMb sdata',siz));
save([pwd '\klustest3\' sdata.combined_name '\' sdata.combined_name '_sdata.mat'],'sdata','-v7.3'); % save session data
save([pwd '\klustest3\' sdata.combined_name '\' sdata.combined_name '_mtint.mat'],'mtint','-v7.3'); % save mtint file
            







