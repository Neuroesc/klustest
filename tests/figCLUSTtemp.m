function figCLUSTtemp(sdata,uci,part_now,fig_vis,save_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION  short descr.
% long descr.
%
% USAGE:
%         [out] = template(in,in2)
%
% INPUT:
%         sdata - sdata structure from klustest
%         uci - unique cell identifier, the function will run on this cell
%         pname - part name to run on
%         fig_vis - (optional), 'on' to show figures, 'off' to hide figures as they are created, default is 'off'
%         save_fig - (optional), 1 to save figures in .fig files, 0 otherwise, default is 0
%
% See also: KLUSTEST

% HISTORY:
% version 1.0.0, Release 31/03/17 Initial release
% version 2.0.0, Release 08/08/18 overhauled for klustest update
% version 2.0.0, Release 08/08/18 added overdispersion and better speed graph

%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
close all
% deal with input variables
inps = {'fig_vis','save_fig'};
vals = {'''on''','0'};
for ff = 1:length(inps)
    if ~exist(inps{ff},'var')
        eval([inps{ff} '=' vals{ff} ';']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
% open/create figure
fig_overall = figure('visible',fig_vis,'Units','pixels','Position',[100, 100, 1610, 800]);
set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white
colormap(jet(256)); % to make sure the colormap is not the horrible default one
fsiz = 6; % the fontsize for different texts
flnw = 0.5; % the line width for different plots

% add an annotation to the figure with some important info
ann_str = sprintf('Cell: %s, Part: %s, Analysed: %s, Auto cell classification: %s',uci,part_now,datestr(now,'yyyy-mm-dd-HH-MM-SS'),sdata.(uci).(part_now).cell_type);
annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',8,'LineStyle','none','interpreter','none');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spike plot with black lines for path and red dots for spikes
axes('Units','pixels','Position',[30,560,200,200])
    % position and spike data
    part_duration = sdata.(part_now).duration;    
    pox = sdata.(part_now).pox;
    poy = sdata.(part_now).poy;
    pot = sdata.(part_now).pot;
    spx = sdata.(uci).(part_now).spx;
    spy = sdata.(uci).(part_now).spy;
    spt = sdata.(uci).(part_now).spt;

    % plot position data, excluding pieces not included in this part
    tvals = sdata.part_config.(part_now).times;
    for tt = 1:length(tvals(:,1))
        tindx = pot > tvals(tt,1) & pot < tvals(tt,2);
        plot(pox(tindx),poy(tindx),'k')
    end
    hold on;

    % plot spikes after so they are all on top
    for tt = 1:length(tvals(:,1))
        sindx = spt > tvals(tt,1) & spt < tvals(tt,2);
        plot(spx(sindx),spy(sindx),'ro','MarkerFaceColor','r','MarkerSize',2)
    end
    daspect([1 1 1])
    axis xy off square tight
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
    title(sprintf('%d spikes (%.2f Hz)',numel(spx),numel(spx)/part_duration),'FontSize',fsiz,'FontWeight','normal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% firing rate map
axrt = axes('Units','pixels','Position',[30,330,200,200]);
    % ratemap and its metrics
    ratemap = double(sdata.(uci).(part_now).ratemap);
    skaggs = sdata.(uci).(part_now).spatial_measures.spatial_information;
    spars = sdata.(uci).(part_now).spatial_measures.sparsity;
    MMF = sdata.(uci).(part_now).spatial_measures.mean_method_focus;
    MI = sdata.(uci).(part_now).spatial_measures.mutual_info;

    % plot the ratemap
    im = imagesc(ratemap);
    set(im,'alphadata',~isnan(ratemap));
    title('Ratemap')
    daspect([1 1 1])
    axis xy off tight square
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);   
    title(sprintf('SI: %.2f, Sp: %.2f, MI: %.2f, MMF: %.2f',skaggs,spars*100,MI,MMF),'FontSize',fsiz,'FontWeight','normal');
    tt = text(-0.07,0.15,sprintf('%c: %.2f (Hz), %c: %.2f (Hz), %c: %.2f (Hz)',char(708),nanmax(ratemap(:)),char(181),nanmean(ratemap(:)),char(709),nanmin(ratemap(:))),'FontSize',fsiz,'Units','normalized');
    set(tt,'rotation',90);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dwell time heatmap
axdt = axes('Units','pixels','Position',[30,100,200,200]);
    % get the dwellmap
    dwellmap = double(sdata.(uci).(part_now).dwellmap);

    % plot the dwellmap
    im = imagesc(dwellmap);
    colormap(axdt,hot(256))    
    set(im,'alphadata',logical(dwellmap));
    daspect([1 1 1])
    axis xy off square tight
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
    title(sprintf('%.2fs (%.2f mins)',part_duration,part_duration/60),'FontSize',fsiz,'FontWeight','normal');  
    tt = text(-0.07,0.15,sprintf('%c: %.2f (s), %c: %.2f (s), %c: %.2f (s)',char(708),nanmax(dwellmap(:)),char(181),nanmean(dwellmap(:)),char(709),nanmin(dwellmap(:))),'FontSize',fsiz,'Units','normalized');
    set(tt,'rotation',90);    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spikes vs time plot (running along bottom of figure)
axst = axes('Units','pixels','Position',[30,35,1550,45]);
    % get the spike time data and bin it
    bstvals = 0:1:sdata.session_duration; % vector of 1s time points at which we should calculate spike probability
    spiketime = sdata.(uci).spike_times;
    [bspikes,~] = histc(spiketime,bstvals);
    
    % plot this as an image
    im = imagesc(bspikes');
    colormap(axst,flipud(gray(256)))
    ax = gca;
    ax.YTick = [];
    ax.TickLength = [0.001, 0.001];
    ylabel('Spikes')
    hold on
    axis xy on
    xlabel('Time (s)')
    tvals = sdata.part_config.(part_now).times;
    for tt = 1:length(tvals(:,1))
        line([tvals(tt,1) tvals(tt,2)],[0.5 0.5],'Color','r','LineWidth',4);
    end
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % set the axes to be on top, set their thickness to flnw, set all fonts to size fsiz and make sure the axis and text are black

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% grid autocorrelation
axgs = axes('Units','pixels','Position',[275,560,200,200]);
    % get the grid metrics
    automap = double(sdata.(uci).(part_now).grid_autocorrelation);    
    grid_score = sdata.(uci).(part_now).grid_score;
    grid_wavelength = sdata.(uci).(part_now).grid_metrics.wavelength;
    grid_orientation = sdata.(uci).(part_now).grid_metrics.orientation;

    % plot the grid autocorrelogram
    msk = sdata.(uci).(part_now).grid_metrics.ring_mask;
    imc = imagesc(automap);
    if ~isnan(msk)
        msk(~msk) = 0.5;
        set(imc,'alphadata',msk);    
    end
    title(sprintf('G: %.2f, Wavelength: %.2f, Orientation: %.2f',grid_score,grid_wavelength,grid_orientation),'FontSize',fsiz,'FontWeight','normal');
    caxis([0 nanmax(automap(:))])
    daspect([1 1 1]);
    axis xy off square tight
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);             
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Place field map    
axpf = axes('Units','pixels','Position',[275,330,200,200]);
    % get the place field data
    fieldd = sdata.(uci).(part_now).field_data;
    overD = sdata.(uci).(part_now).over_dispersion;

    % create the binarized ratemap
    bmap = single(imbinarize(ratemap,max([sdata.config.minfr, sdata.config.frcut*nanmax(ratemap(:))]) ));
    bmap(bmap>0) = .75;
    bmap(isnan(ratemap)) = NaN;
    
    % plot the fields and their areas
    im = imagesc(bmap);
    set(im,'alphadata',~isnan(bmap));
    daspect([1 1 1])
    axis xy off tight square
    title(sprintf('%d fields, %carea: %.1fpx, %cSNR: %.1f',length(fieldd.Area(:)),char(181),nanmean(fieldd.Area(:)),char(181),nanmean(fieldd.signal_to_noise(:))),'FontSize',fsiz,'FontWeight','normal');
    colormap(axpf,bone(128))
    caxis([0 1]);
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    text(fieldd.WeightedCentroid(:,1),fieldd.WeightedCentroid(:,2),cellstr(num2str(fieldd.Area(:))),'Color','r','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fsiz);
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interspike intervals
axes('Units','pixels','Position',[275,115,200,185])
    yyaxis left
    adist = sdata.(uci).(part_now).isi_data.adist; % actual isi histogram   
    if ~all(isnan(adist))
        bar(adist(:,1),adist(:,2),1,'k');
        hold on
        bar(-adist(:,1),adist(:,2),1,'k');
        ax = gca;
        ylabel('Spikes')
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); 

        yyaxis right
        fdist = sdata.(uci).(part_now).isi_data.fdist; % fitted isi histogram   
        fx = fdist(:,1);
        plot(fdist(:,1),fdist(:,2),'r','LineStyle','-');
        hold on
        plot(-fdist(:,1),fdist(:,2),'r','LineStyle','-');

        ax = gca;
        ax.YTick = [];
        hmax = sdata.(uci).(part_now).isi_data.half_max;
        ps = sdata.(uci).(part_now).isi_data.hwidth_ps;
        line([fx(ps(2)) fx(ps(3))],[hmax hmax],'Color','r','LineWidth',1.5)
        line([fx(ps(1)) fx(ps(1))],ax.YLim,'Color','r','LineWidth',1)

        box on
        xlabel('Interspike interval (ms)') % label x axis
        astr = sprintf('fwhm: %.2f, sd: %.2f',sdata.(uci).(part_now).isi_data.fwhmx,sdata.(uci).(part_now).isi_data.stdev);
        title(astr,'FontSize',fsiz,'FontWeight','normal');
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); 
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phase preference plot     
axphe = axes('Units','pixels','Position',[520,665,275,95]);
    % get the phase preference data
    phprobs = sdata.(uci).(part_now).spike_phase_ksdensity(:,1);
    swav = sdata.(uci).(part_now).spike_phase_ideal;
    yi = sdata.(uci).(part_now).spike_phase_binned;
    aip = (deg2rad(-180:5:540))'; % doing this means we have bins symmetrical around zero

    % plot the kernel smoothed density estimate in red and the theta phase in blue dotted
    yyaxis left
    plot(aip,phprobs,'r','LineWidth',flnw)
    ax = gca;
    ax.YTick = [];
    ylabel('Probability')
    
    yyaxis right
    plot(aip,swav,'b','LineWidth',flnw);    
    ax = gca;
    ax.YTick = [];
    
    ax = gca;
    ax.YAxis(1).Limits = [min(phprobs) max(phprobs)];
    ax.YAxis(2).Limits = [min(swav) max(swav)];    
    ax.XAxis.Limits = [min(aip) max(aip)];
    text(0.1,1.1,'Theta','Color','b','FontSize',fsiz,'Units','normalized');
    text(0.8,1.1,'KDE','Color','r','FontSize',fsiz,'Units','normalized');    
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  

axphe = axes('Units','pixels','Position',[520,560,275,100]);
    % plot the firing rate by phase
    bar(aip,yi,1,'k');
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  
    ax = gca;
    ax.YColor = 'k';
    ylabel('Spikes')
    xlabel(sprintf('Theta Phase (%c)',char(176)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spike autocorrelation - theta plot           
axes('Units','pixels','Position',[520,340,275,170])
    % get the theta data
    theta_data = sdata.(uci).(part_now).theta_train_long;
    corrdata = theta_data(:,1);
    tms = theta_data(:,2);
    corrfilt = theta_data(:,3);
    theta_index = sdata.(uci).(part_now).theta_index;
    theta_ratio = sdata.(uci).(part_now).theta_ratio;
    theta_skip = sdata.(uci).(part_now).theta_skip;

    yyaxis right
    plot(tms,corrfilt,'Color','r','linewidth',1);
    ax = gca;
    ax.YTick = [];
    ax.YColor = 'k';

    yyaxis left
    bar(tms,corrdata,1,'k');
    ax = gca;
    ax.YColor = 'k';
    axis([tms(1) tms(end) ax.YLim]);
    xlabel('Time lag (ms)') % label x axis
    ylabel('Probability') % label y axis
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
    astr = sprintf('%c index: %.2f, %c ratio: %.2f, %c skip: %.2f',char(952),theta_index,char(952),theta_ratio,char(952),theta_skip);
    title(astr,'FontSize',fsiz,'FontWeight','normal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spike autocorrelation - refractory period plot   
axes('Units','pixels','Position',[520,120,275,170])
    % get the refractory period info
    burst_index = sdata.(uci).(part_now).burst_index;
    burst_length_mean = sdata.(uci).(part_now).burst_length_mean;
    tms2 = sdata.(uci).(part_now).refractory_period(:,2);
    corrdata2 = sdata.(uci).(part_now).refractory_period(:,1);
    tau_r = sdata.tau_r;

    bar(tms2,corrdata2,0.9,'k');
    ax = gca;
    axis([tms2(1) tms2(end) ax.YLim]);
    xlabel('Time lag (ms)') % label x axis
    ylabel('Probability')
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);     

    astr = sprintf('burst index: %.2f, burst %clength: %.2f, RPVs: %d, RPVp: %.1f, fpRate: %.2f',burst_index,char(956),burst_length_mean,sdata.(uci).(part_now).rpv_total,sdata.(uci).(part_now).rpv_proportion,sdata.(uci).(part_now).rpv_false_positive1);
    title(astr,'FontSize',fsiz,'FontWeight','normal');

    hold on
    plot([-tau_r; -tau_r],ax.YLim,'r','LineWidth',1)
    plot([tau_r; tau_r],ax.YLim,'r','LineWidth',1) 
    plot([-6; -6],ax.YLim,'r:','LineWidth',0.5)
    plot([6; 6],ax.YLim,'r:','LineWidth',0.5) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Waveform plot
axes('Units','pixels','Position',[840,560,200,200])
    % get the wave data
    wavtime = -200:20:780;
    maxwavs = sdata.(uci).(part_now).waveform_max;
    [ax_high,mval] = nanmax(maxwavs);
    ax_high = ax_high + max(max(sdata.(uci).(part_now).waveform_stdv{mval}));
    ax_low = sdata.(uci).(part_now).waveform_min(mval) - max(max(sdata.(uci).(part_now).waveform_stdv{mval}));
    
    % plot the mean and SD of the waveform with the highest amplitude
    [hl,hp] = boundedline(wavtime,sdata.(uci).(part_now).waveform_mean{mval},sdata.(uci).(part_now).waveform_stdv{mval},'-k');
    set(hl,'Color','r') % line color
    set(hl,'LineStyle','-') % line style
    set(hl,'LineWidth',1) % line width
    set(hp,'FaceColor','b') % color of area
    set(hp,'FaceAlpha',0.5) % transparency of area
    ax = gca;
    ax.XLim = [-200 780];
    ax.YLim = [ax_low ax_high];
    box on
    grid on
    axis on square tight xy
    xlabel('Time (ms)')
    ylabel(sprintf('Amplitude (%cV)',char(956)))
    
    astr = sprintf('Peak: %.1f%cV, Width: %.fms, SNR: %.2f',sdata.(uci).(part_now).waveform_max(mval),char(956),sdata.(uci).(part_now).waveform_width(mval),sdata.(uci).(part_now).channel_snr(mval));
    title(astr,'FontSize',fsiz,'FontWeight','normal');
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);     

    % plot all 4 waveforms in smaller plots
    wpos = [1045,710,50,50; 1045,660,50,50; 1045,610,50,50; 1045,560,50,50];
    for ww = 1:4
        axes('Units','pixels','Position',wpos(ww,:))
        wavtime = -200:20:780;

        if ~all(sdata.(uci).(part_now).waveform_mean{ww}) % if the waveform is just zeros
            ax = gca;
            ax.XLim = [-200 780]; 
            ax.YLim = [ax_low ax_high];    
            ax.XTick = [];
            ax.YTick = [];
            box on   
            axis square        
            text(-100,0,'Grounded','FontSize',7)
            continue
        end % if ~all(sdata.(uci).(part_now).waveform_mean{ww})

        [hl,hp] = boundedline(wavtime,sdata.(uci).(part_now).waveform_mean{ww},sdata.(uci).(part_now).waveform_stdv{ww},'-k');

        % we want the maximum waveform to be different
        if ww == mval
            set(hl,'Color','r') % line color
            set(hl,'LineStyle','-') % line style
            set(hl,'LineWidth',1) % line width
            set(hp,'FaceColor','b') % color of area
            set(hp,'FaceAlpha',.5) % transparency of area        
        else
            set(hl,'Color','w') % line color
            set(hl,'LineStyle','-') % line style
            set(hl,'LineWidth',1) % line width
            set(hp,'FaceColor','k') % color of area
            set(hp,'FaceAlpha',.5) % transparency of area
        end % if ww == mval

        ax = gca;
        ax.XLim = [-200 780]; 
        ax.YLim = [ax_low ax_high];    
        ax.XTick = [];
        ax.YTick = [];
        box on   
        grid on
        axis square
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);             
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mahalanobis distance, cluster quality plot
axes('Units','pixels','Position',[840,340,255,170]);
    % determine cluster and tetrode
    sst = strsplit(uci,'_');
    tet = str2double(sst{3}(regexp(sst{3},'\d')));
    tetstr = ['t' num2str(tet)];
    clu = str2double(sst{4}(regexp(sst{4},'\d')));

    % cluster space data
    ndat = single(sdata.(tetstr).fetdata.noise_dists{clu});
    cdat = single(sdata.(tetstr).fetdata.clust_dists{clu});
    isod = sdata.(tetstr).fetdata.isolation_distances(clu);
    lratio = sdata.(tetstr).fetdata.lratios(clu);

    % plot this info
    if ~isempty(cdat) && ~isempty(ndat)
        % get maximum x value
        max_d = ceil(max([max(ndat);max(cdat)]))+1;
        if (isnan(max_d)) % to cope with missing channel data
            max_d = 1;
        end 

        xi = linspace(0,max_d,10000); % vector of values where we want to estimate ksdensity
        [vals1,~,~] = ksdensity(cdat,xi,'Support',[-1 max_d],'Kernel','epanechnikov','Function','pdf','Bandwidth',0.05);
        [vals2,~,~] = ksdensity(ndat,xi,'Support',[-1 max_d],'Kernel','epanechnikov','Function','pdf','Bandwidth',0.05);

        a1 = area(xi,vals1);
        set(a1,'FaceColor','b')
        alpha(.25)
        hold on
        a2 = area(xi,vals2);
        set(a2,'FaceColor','k')
        alpha(.25)

        set(gca,'XScale','log');
        set(gca,'Xlim',[0,max_d]);
        xlabel('Log(Distance)') % label x axis
        ylabel('Probability') % label y axis
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);  
        astr = sprintf('IsoD: %.2f, Lratio: %.2f',isod,lratio);
        title(astr,'FontSize',fsiz,'FontWeight','normal');

    else
        text(0.1,0.5,'Not computed');
    end 
    box on
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster space plot  
axes('Units','pixels','Position',[840,120,255,170])
    % cluster and cluster space info
    fetdata = sdata.(tetstr).fetdata;
    fdata = fetdata.fetdata;
    nfets = fetdata.nfeatures;
    clusters_now = unique(fetdata.cluster_labels);
    clus_count = numel(clusters_now);
    clus_cut = fetdata.cluster_labels;
    tspikes = fetdata.nspikes;

    % find the 2 channels with the maximum waveforms
    [~,dindx] = sort(maxwavs,2,'descend'); % find the order of these values
    mch1 = dindx(1); % the channel with the maximum amplitude
    mch2 = dindx(2); % the channel with the second maximum amplitude                    
    start_cols = 1:nfets:nfets*4; % the starting column for each channel

    % get the feature data for the two channels
    d1 = fdata(:,start_cols(mch1)); % get this feature data for this channel
    d2 = fdata(:,start_cols(mch2)); % get this feature data for this channel

    % plot the noise in grey and then the cell in red
    colplot = [0.5 0.5 0.5 0.5];
    plot(d1,d2,'.','MarkerSize',3,'color',colplot); % plot the clusters
    hold on
    cindx = find(clus_cut == clu);
    plot(d1(cindx),d2(cindx),'.','MarkerSize',3,'color','r'); % re-plot the current cluster in red to make sure it is on top and stands out

    axis xy on
    astr = sprintf('Ch%d vs Ch%d: %d clusters, %d spikes (%.1f%% of total %d)',mch1,mch2,clus_count,numel(spx),numel(spx)/tspikes*100,tspikes);
    title(astr,'FontSize',fsiz,'FontWeight','normal');
    xlabel(sprintf('PC1 Ch%d',mch1))
    ylabel(sprintf('PC1 Ch%d',mch2))
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);      
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% head direction linear plot
axhd2 = axes('Units','pixels','Position',[1140,340,200,170]);   
    % HD data
    hd_type = sdata.config.hd_type;
    if strcmp(hd_type,'histogram')
        hd1 = sdata.(uci).(part_now).hd_session;
        hd3 = sdata.(uci).(part_now).hd_cell;
        hd_c = sdata.(uci).(part_now).hd_frate;
        rayleigh = sdata.(uci).(part_now).hd_rayleigh;
        mx2 = sdata.(uci).(part_now).hd_maximum;
        mn2 = sdata.(uci).(part_now).hd_mean; % add data to structure                
        sd2 = sdata.(uci).(part_now).hd_stdev; % add data to structure                
    else
        hd1 = sdata.(uci).(part_now).hd_density_session;
        hd3 = sdata.(uci).(part_now).hd_density_cell;
        hd_c = sdata.(uci).(part_now).hd_density_frate;
        rayleigh = sdata.(uci).(part_now).hd_density_rayleigh;
        mx2 = sdata.(uci).(part_now).hd_density_maximum;
        mn2 = sdata.(uci).(part_now).hd_density_mean; % add data to structure                
        sd2 = sdata.(uci).(part_now).hd_density_stdev; % add data to structure                
    end
    ai = linspace(0,2*pi,sdata.config.hd_bins)'; % angles for binning   
    
    % plot the data as a linear (side on) graph
    a1 = area(rad2deg(ai),hd1);
    set(a1,'FaceColor','k')
    alpha(.25)    
    hold on
    a2 = area(rad2deg(ai),hd3);
    set(a2,'FaceColor','b')
    alpha(.25)        
    
    ax = gca;
    ax.XLim = [0 360];
    ax.XTick = [0:60:360];
    xlabel(sprintf('HD (%c)',char(176)))
    if strcmp(hd_type,'histogram')
        ylabel('Firing Rate (Hz)')        
    else
        ylabel('Normalized probability')
    end
    title(sprintf('r: %.2f, %c: %.2f, %c: %.2f, %s: %.2f',rayleigh,char(708),mx2,char(956),mn2,char(963),sd2),'FontSize',fsiz,'FontWeight','normal');
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);      
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% firing rate by speed
axvs = axes('Units','pixels','Position',[1140,120,200,170]);
    % speed data
    speed_time = sdata.(uci).(part_now).speed_time;
    sscore = sdata.(uci).(part_now).speed_score;
    sslope = sdata.(uci).(part_now).speed_slope;
    sintpt = sdata.(uci).(part_now).speed_intercept;
    scurve = sdata.(uci).(part_now).speed_curve;

    % plot this relationship
    yyaxis left
    plot(scurve(:,1),scurve(:,2),'k');
    hold on;
    scatter(scurve(:,1),scurve(:,2),20,'k')
    hold on
    rr = refline(sslope,sintpt);
    set(rr,'Color','r');
    title(sprintf('Sscore: %.2f, slope: %.2f, intercept: %.2f',sscore,sslope,sintpt),'FontSize',fsiz,'FontWeight','normal');
    ax = gca;
    ax.YLim = [0 max(scurve(1:25,2))+0.1];
    xlabel('Speed (cm/s)')
    ylabel('Firing Rate (Hz)')
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k');

    yyaxis right
    plot(scurve(:,1),speed_time(:),'b');
    hold on
    scatter(scurve(:,1),speed_time(:),20,'b')
    ylabel('Time (s)')
    ax = gca;
    ax.YLim = [0 max(speed_time(1:25,1))+10];
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','b');

    axis on tight xy
    ax.XLim = [0 50];  
    box on
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% overdispersion z distribution
axov = axes('Units','pixels','Position',[1390,560,200,200]);
    % get the OD z distribution
    zdist = sdata.(uci).(part_now).over_dispersion_z;

    % get the KDE for this
    edg = -10:0.1:10;
    fz = ksdensity(zdist,edg,'kernel','epanechnikov','Function','pdf');
    fg = gaussmf(edg,[1,0]);
    fg = fg .* nanmax(fz(:));

    plot(edg,fz,'k')
    hold on
    plot(edg,fg,'r')

    ax = gca;
    ax.XLim = [-10 10];  
    xlabel('Standardized rate')
    ylabel('Probability')
    box on
    title(sprintf('%c: %.1f, %c{^2}: %.1f, N: %d ',956,nanmean(zdist),963,overD,numel(zdist)),'FontSize',fsiz,'FontWeight','normal');
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cell classification text
axtxt = axes('Units','pixels','Position',[1390,120,200,170]);
    % get the cell's firing statistics
    frate = sdata.(uci).(part_now).frate;
    [wpeak,pindx] = nanmax(sdata.(uci).(part_now).waveform_max);
    wwide = sdata.(uci).(part_now).waveform_width(pindx);
    skagg = sdata.(uci).(part_now).spatial_measures.spatial_information;
    gscore = sdata.(uci).(part_now).grid_score;
    rvect = sdata.(uci).(part_now).hd_rayleigh;

    % place text in axis
    astr = sprintf('Mean firing rate: %.1f Hz\nWaveform width: %.1f ms\nSkaggs spatial information: %.1f b/s\nGrid score: %.1f\nRayleigh vector: %.1f\nAuto classified as: %s',frate,wwide,skagg,gscore,rvect,sdata.(uci).(part_now).cell_type);
    text(0.1,0.5,astr,'Units','normalized','FontSize',8,'FontWeight','normal')
    axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% head direction polar plot
% WARNING: for some reason this plot must be done last, something in the mmpolar code messes up the axis calling
axhd = axes('Units','pixels','Position',[1140,560,200,200]);
    % polar plot of HD data
    set(fig_overall,'CurrentAxes',axhd);
    mmp = mmpolar(ai,hd1,'k',ai,hd3,'b','FontSize',fsiz,'Grid','on','RGridVisible','off','RTickVisible','off','TTickDelta',20,'RTickLabelVisible','on','TTickLabelVisible','on');
    set(mmp(1),'LineWidth',0.5);
    set(mmp(2),'LineWidth',0.5);
    p1 = patch(get(mmp(1),'XData'),get(mmp(1),'YData'),'k','FaceAlpha',0.2);
    p2 = patch(get(mmp(2),'XData'),get(mmp(2),'YData'),'b','FaceAlpha',0.5);
    hl = legend([p1,p2],{'Session','Cell'},'Units','normalized','Position',[.655,.655,.15,.05]);
    legend boxoff
    axis square
    set(hl,'FontSize',6);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the figure    
figfile = [pwd '\klustest\' sdata.combined_name '\figures\'];
[~,~,~] = mkdir(figfile);
print(fig_overall,'-dpng','-r150',[figfile uci '_' part_now '.png'])
if save_fig
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')');
    savefig(fig_overall,[figfile uci '_' part_now '.fig'],'compact');
end
close(fig_overall);                                          

















































