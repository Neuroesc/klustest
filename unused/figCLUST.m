function figCLUST(mtint,pdata,sdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figCLUST  a basic figure function
% This function just stores the plotting for klustest. This figure gives all the information about a cluster
% in only one part (i.e. waveforms, ratemap, refractory period)
%
% USAGE:
%         figCLUST(sdata,uci,part_now,fig_vis,save_fig)
%
% INPUT:
%         sdata - sdata structure from klustest
%         uci - unique cell identifier, the function will run on this cell
%         part_now - part name to run on
%         fig_vis - (optional), 'on' to show figures, 'off' to hide figures as they are created, default is 'off'
%         save_fig - (optional), 1 to save figures in .fig files, 0 otherwise, default is 0
%
% See also: KLUSTEST

% HISTORY:
% version 1.0.0, Release 31/03/17 Initial release
% version 2.0.0, Release 08/08/18 overhauled for klustest update
% version 2.0.0, Release 08/08/18 added overdispersion and better speed graph
% version 3.0.0, Release 10/04/19 overhauled for klustest update
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

srow = size(sdata,1);
if ~sdata.spikes(srow) || isempty(sdata.spikes(srow)) || isnan(sdata.spikes(srow))
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
    fig_clust = figure('visible',fig_vis,'Units','pixels','Position',[50, 50, 1610, 920]); % open/create figure
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    fsiz = 6; % the fontsize for different texts
    flnw = 0.5; % the line width for different plots

    % add an annotation to the figure with some important info
    part_now = sdata.part{srow};
    ann_str = sprintf('Cell: %d, Rat: %d, Date: %d, Tetrode: %d, Cluster: %d, Part: %s, Analysed: %s, Auto cell classification: %s',sdata.uci(srow),sdata.rat(srow),sdata.date(srow),sdata.tet(srow),sdata.clu(srow),part_now,datestr(now,'yyyy-mm-dd-HH-MM-SS'),sdata.cell_type{srow});
    annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',8,'LineStyle','none','interpreter','none');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Position and spike plot (black lines & red dots)
    ax_ps = axes('Units','pixels','Position',[30,560,200,200]);
        % position and spike data (already cut to this part)
        part_duration = sdata.duration(srow);  
        ppox = pdata.(part_now).pox;
        ppoy = pdata.(part_now).poy;
        ppot = pdata.(part_now).pot; 
        pspx = pdata.(part_now).pox(sdata.spike_index{srow});
        pspy = pdata.(part_now).poy(sdata.spike_index{srow});
        pspt = sdata.spike_times{srow};
        
        % plot position data, excluding pieces not included in this part
        % by inserting NaNs between intervals we can plot this as one line, which saves on memory
        dindax = abs([0; diff(ppot)])>0.1;
        pos_plot = [ppox ppoy];
        pos_plot(dindax,:) = NaN;
        plot(pos_plot(:,1),pos_plot(:,2),'Color',[.5 .5 .5]); hold on;

        % plot spikes after position data so they are all on top
        plot(pspx,pspy,'ro','MarkerFaceColor','r','MarkerSize',3)    
        
        % additional settings
        daspect([1 1 1])
        axis xy off square tight
        text(.5,1.05,sprintf('%d spikes (%.2f Hz)',numel(pspx),numel(pspx)/part_duration),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        td = nansum(sqrt(sum([diff(ppox).^2,diff(ppoy).^2],2)));
        tt = text(-0.07,0.5,sprintf('Path dist: %.2fm, Width: %.2fm, Height: %.2fm',td/100,range(ppox(:))/100,range(ppoy(:))/100),'FontSize',fsiz+1,'Units','normalized','HorizontalAlignment','center');
        set(tt,'rotation',90);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% 15 spike plots showing consistency over time
    % work out how much time/samples must be shown per plot
    total_samples = numel(ppot);
    samples_per_plot = max([60/(1/50) ceil(total_samples / 15)]); % a minimum of 60s per plot
    time_index = 1:samples_per_plot:total_samples;

    axy = 786;
    axx = [30 140 250 360 470 580 690 800 910 1020 1130 1240 1350 1460];
    axs = [100 100];
    for tt = 1:(length(time_index)-1)
        ax_ps1 = axes('Units','pixels','Position',[axx(tt),axy,axs]); 
            plot(pos_plot(time_index(tt):time_index(tt+1),1),pos_plot(time_index(tt):time_index(tt+1),2),'Color',[.5 .5 .5]); hold on;
            sindx = pspt > ppot(time_index(tt)) & pspt < ppot(time_index(tt+1));
            plot(pspx(sindx),pspy(sindx),'ro','MarkerFaceColor','r','MarkerSize',2)    

            ax = gca;
            ax.XLim = ax_ps.XLim;
            ax.YLim = ax_ps.YLim;            
            ax.YTick = [];
            ax.XTick = [];
            set(gca,'LineWidth',flnw,'layer','top','XColor',[.9 .9 .9],'YColor',[.9 .9 .9]);           
            text(.5,1.05,sprintf('%.1f - %.1fs',ppot(time_index(tt)),ppot(time_index(tt+1))),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Spike rasterplot running along bottom of figure
    axsr = axes('Units','pixels','Position',[30,35,1550,45]);
        % get the spike time data and bin it
        bsize = 5;
        edg = 0:bsize:pdata.duration; % vector of 1s time points at which we should calculate spike probability
        spt = mtint.tetrode(sdata.tet(srow)).ts( mtint.tetrode(sdata.tet(srow)).cut==sdata.clu(srow) );
        f = histcounts(spt,edg);

        % plot this as an image
        im = imagesc(f);
        colormap(axsr,flipud(gray(256)))
        ax = gca;
        ax.YTick = [];
        ax.TickLength = [0.001, 0.001];
        ylabel('Spikes')
        hold on
        axis xy on
        xlabel('Time (s)')
        tvals = pdata.(part_now).times;
        
        for tt = 1:length(tvals(:,1))
            line([knnsearch(edg(:),tvals(tt,1)) knnsearch(edg(:),tvals(tt,2))],[0.5 0.5],'Color','r','LineWidth',4);
        end
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % set the axes to be on top, set their thickness to flnw, set all fonts to size fsiz and make sure the axis and text are black 
        set(ax,'xticklabel',num2str(get(gca,'xtick')'.*bsize,'%.f'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Firing rate map
    axrt = axes('Units','pixels','Position',[30,330,200,200]);
        ratemap = sdata.ratemap{srow};
        spati = sdata.spatial_info_bsec(srow);
        spars = sdata.sparsity(srow);
        coher = sdata.spatial_coherence(srow);

        im = imagesc(ratemap);
        set(im,'alphadata',~isnan(ratemap));
        daspect([1 1 1])
        caxis([0 nanmax(ratemap(:))])        
        axis xy off tight square
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);   
        
        text(.5,1.05,sprintf('SI: %.2f, Sp: %.2f, Cohe: %.2f',spati,spars,coher),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        tt = text(-0.07,0.5,sprintf('Max: %.2f, Mean: %.2f, Median: %.2f',nanmax(ratemap(:)),nanmean(ratemap(:)),nanmedian(ratemap(:))),'FontSize',fsiz+1,'Units','normalized','HorizontalAlignment','center');
        set(tt,'rotation',90);
        
        axp = get(gca,'Position');
        cc = colorbar;
        set(gca,'Position',axp);
        set(cc,'Position',get(cc,'Position')+[0 0 -0.004 0])
        title(cc,'Hz','FontSize',fsiz)
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Dwell time map
    axdt = axes('Units','pixels','Position',[30,100,200,200]);
        dwellmap = sdata.dwellmap{srow};

        im = imagesc(dwellmap);
        colormap(axdt,hot(256))    
        set(im,'alphadata',logical(dwellmap));
        daspect([1 1 1])
        caxis([0 nanmax(dwellmap(:))])
        axis xy off square tight
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
        title(sprintf('%.2fs (%.2f mins)',part_duration,part_duration/60),'FontSize',fsiz,'FontWeight','normal');  
        tt = text(-0.07,0.5,sprintf('Max: %.2f, Mean: %.2f, Median: %.2f',nanmax(dwellmap(:)),nanmean(dwellmap(:)),nanmedian(dwellmap(:))),'FontSize',fsiz+1,'Units','normalized','HorizontalAlignment','center');
        set(tt,'rotation',90);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Grid autocorrelation
    axgs = axes('Units','pixels','Position',[275,560,200,200]);
        amap = sdata.grid_autocorrelation{srow};    
        msk = sdata.grid_mask{srow};
        
        imc = imagesc(amap);
        set(imc,'alphadata',double(msk)+0.3);    
        daspect([1 1 1])
        caxis([0 1])
        axis xy off square tight
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
        text(.5,1.05,sprintf('G: %.2f, W: %.2f, O: %.2f',sdata.grid_score(srow),sdata.grid_wavelength(srow),sdata.grid_orientation(srow)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Place fields
    axpf = axes('Units','pixels','Position',[275,330,200,200]);
        lowfr = max([pdata.config.minfr, pdata.config.frcut*nanmax(ratemap(:))]);
        highfr = max([pdata.config.minfr, nanmax(ratemap(:))+0.01]);
        
        im = imshow(ratemap,[lowfr highfr]); hold on;
        set(im,'AlphaData',~isnan(ratemap))
        daspect([1 1 1])
        axis xy off square tight        
        
        for ff = 1:sdata.place_fields(srow)
            cnt = sdata.field_centroids{srow}(ff,:);
            a = sdata.field_maj_lengths{srow}(ff,:);
            b = sdata.field_min_lengths{srow}(ff,:);
            o = deg2rad(sdata.field_orientations{srow}(ff,:));

            % plot ellipse
            t = linspace(0,2*pi,100);
            x = cnt(1) + a*cos(t)*cos(o) - b*sin(t)*sin(o);
            y = cnt(2) + b*sin(t)*cos(o) + a*cos(t)*sin(o);
            plot(x,y,'Color',autumn(ff))
            text(cnt(1),cnt(2),num2str(ff),'Color',autumn(ff),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
        end
        
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
        text(.5,1.05,sprintf('Fields: %d',sdata.place_fields(srow)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Interspike interval       
    axes('Units','pixels','Position',[275,130,200,170])
        xi = 0:0.5:50;
        idist = sdata.isi_dist{srow};
        fdist = sdata.isi_fdist{srow};
        fdist = fdist./nanmax(fdist(:)) .* nanmax(idist); 
        
        bar(xi,idist,1,'k'); hold on;
        area(-xi,fdist,'FaceColor','k')
       
        ax = gca;
        ax.XLim = [-50 50];
        hmax = sdata.isi_half_width(srow,1);
        ps = sdata.isi_half_width(srow,2:end);
        line(-[ps(2) ps(3)],[hmax hmax]./nanmax(sdata.isi_fdist{srow}(:)).*nanmax(idist),'Color','r','LineWidth',1.5)
        line(-[ps(1) ps(1)],ax.YLim,'Color','r','LineWidth',1)
        
        ylabel('Spikes')
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); 
        xlabel('Interspike interval (ms)') % label x axis
        text(.5,1.05,sprintf('fwhmx: %.1fms, 10ms proportion: %.2f',sdata.isi_fwhmx(srow),nansum(idist(xi<=10))./numel(pspt)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        grid on;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Theta phase preference plot            
    axphe = axes('Units','pixels','Position',[520,590,275,170]);
        phase_dist = sdata.theta_phase_dist{srow};
        ai = reshape(deg2rad(-180:5:540),[],1); % doing this means we have bins symmetrical around zero
        xi = histcents(ai);
        
        bar(xi,phase_dist,1,'k'); hold on;
        ax = gca;
        y1 = ax.YLim;
        si = cos(xi) ./ max(sin(xi)) .* (y1(2)/2) + (y1(2)/2);
        plot(xi,si,'Color','r','LineWidth',1);
        ax.YLim = y1;

        ylabel('Spikes')
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); 
        ax.XTick = [-pi 0 pi 2*pi 3*pi];
        ax.XTickLabel = {'-\pi','0','\pi','2\pi','3\pi'};
        xlabel('Theta phase (rad)') % label x axis
        text(.5,1.05,sprintf('mean: %.1f, r: %.2f',sdata.theta_phase_mean(srow),sdata.theta_phase_r(srow)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Spike autocorrelation 1 - 500ms theta
    axes('Units','pixels','Position',[520,360,275,170])
        autoc = sdata.t500_spike_autocorr{srow};
        tlag = -495:10:495;
        autof = sdata.t500_spike_autocorr{srow};
      
        bar(tlag,autoc,1,'k'); hold on;
        plot(tlag,autof,'Color','r','LineWidth',1)
        
        ax = gca;
        ax.XLim = [-500 500];
        xlabel('Time lag (ms)') % label x axis
        ylabel('Probability') % label y axis
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
        tf = sdata.intrinsic_theta_frequency(srow);
        y1 = ax.YLim;
        la = line([1/tf*1000 1/tf*2000 1/tf*3000 1/tf*4000;1/tf*1000 1/tf*2000 1/tf*3000 1/tf*4000],[ax.YLim',ax.YLim',ax.YLim',ax.YLim'],'Color',[.5 .5 .5]);
        lb = line(-[1/tf*1000 1/tf*2000 1/tf*3000 1/tf*4000;1/tf*1000 1/tf*2000 1/tf*3000 1/tf*4000],[ax.YLim',ax.YLim',ax.YLim',ax.YLim'],'Color',[.5 .5 .5]);   
        uistack(la,'bottom'); uistack(lb,'bottom');
        ax.YLim = y1;
        text(.5,1.05,sprintf('Theta index: %.1f, frequency: %.2fHz, fit: %.2f',sdata.intrinsic_theta_index(srow),sdata.intrinsic_theta_frequency(srow),sdata.intrinsic_theta_fit(srow)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Spike autocorrelation 1 - 25ms refractory period    
    axes('Units','pixels','Position',[520,130,275,170])
        autoc = sdata.t25_spike_autocorr{srow};
        tlag = -24.75:0.5:24.75;
      
        bar(tlag,autoc,1,'k'); hold on;
        ax = gca;
        tau_r = 2; % length of refractory period in ms
        plot([-tau_r; -tau_r],ax.YLim,'r','LineWidth',1)
        plot([tau_r; tau_r],ax.YLim,'r','LineWidth',1) 
        plot([-6; -6],ax.YLim,'r:','LineWidth',0.5)
        plot([6; 6],ax.YLim,'r:','LineWidth',0.5) 
        
        ax = gca;
        ax.XLim = [-25 25];
        xlabel('Time lag (ms)') % label x axis
        ylabel('Probability') % label y axis
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
        text(.5,1.05,sprintf('Burst index: %.1f, RPV: %d (%.2f)',sdata.burst_index(srow),sdata.rpv_total(srow),sdata.rpv_proportion(srow)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.3f'))
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Waveform plots      
    axes('Units','pixels','Position',[840,590,200,170])
        wavtime = -200:20:780;
        maxwavs = sdata.waveform_max(srow,:);
        [maxval,maxidx] = nanmax(maxwavs);
        mwaves = reshape(cell2mat(sdata.waveform_mean(srow,:)),1,50,4);
        swaves = reshape(cell2mat(sdata.waveform_stdv(srow,:)),1,50,4);
        ax_lim = [floor(nanmin(mwaves(:)-swaves(:))./10).*10 ceil(nanmax(mwaves(:)+swaves(:))./10).*10];
        
        % plot the mean and SD of the waveform with the highest amplitude
        [hl,hp] = boundedline(wavtime,mwaves(1,:,maxidx),swaves(1,:,maxidx),'-k');
        set(hl,'Color','r','LineStyle','-','LineWidth',1); set(hp,'FaceColor','b','FaceAlpha',.5);     
        ax = gca;
        ax.XLim = [-200 780];
        ax.YLim = ax_lim;
        box on
        grid on
        xlabel('Time (ms)')
        ylabel(sprintf('Amplitude (%cV)',char(956)))
        text(.5,1.05,sprintf('Peak: %.1f%cV, Width: %.fms, SNR: %.1f',maxval,char(956),sdata.waveform_width(srow,maxidx),sdata.waveform_snr(srow,maxidx)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);     

        % plot all 4 waveforms in smaller plots
        wpos = [1045,710,50,50; 1045,660,50,50; 1045,610,50,50; 1045,560,50,50];
        for ww = 1:4
            axes('Units','pixels','Position',wpos(ww,:))
                [hl,hp] = boundedline(wavtime,mwaves(1,:,ww),swaves(1,:,ww),'-k');

                % we want the maximum waveform to be different
                if ww == maxidx
                    set(hl,'Color','r','LineStyle','-','LineWidth',1); set(hp,'FaceColor','b','FaceAlpha',.5);     
                else
                    set(hl,'Color','k','LineStyle','-','LineWidth',1); set(hp,'FaceColor','k','FaceAlpha',.5);     
                end

                ax = gca;
                ax.XLim = [-200 780]; 
                ax.YLim = ax_lim;
                ax.XTick = [];
                ax.YTick = [];
                box on   
                grid on
                set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);             
        end        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Cluster quality
    axes('Units','pixels','Position',[840,360,255,170]);
        ndat = single(mtint.clu_quality(sdata.tet(srow)).noise_dists{sdata.clu(srow)});
        cdat = single(mtint.clu_quality(sdata.tet(srow)).clust_dists{sdata.clu(srow)});
        isod = sdata.isod(srow);
        lrat = sdata.lratio(srow);

        if ~isempty(cdat) && ~isempty(ndat)
            max_d = nanmax([ceil(max([max(ndat);max(cdat)]))+1, 1]); % maximum x value

            xi = linspace(0,max_d,1000); % vector of values where we want to estimate ksdensity
            [vals1,~,~] = ksdensity(cdat,xi,'Support',[-1 max_d],'Kernel','epanechnikov','Function','pdf');
            [vals2,~,~] = ksdensity(ndat,xi,'Support',[-1 max_d],'Kernel','epanechnikov','Function','pdf');

            a1 = area(xi,vals1,'FaceColor','b','FaceAlpha',.25); hold on;
            a2 = area(xi,vals2,'FaceColor','k','FaceAlpha',.25);

            set(gca,'XScale','log');
            set(gca,'Xlim',[0,max_d]);
            xlabel('Log(Distance)') % label x axis
            ylabel('Probability') % label y axis
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);  
            set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))            
            text(.5,1.05,sprintf('Iso-D: %.2f, L-ratio: %.2f',isod,lrat),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
            legend([a1 a2],{'Cluster','Noise'},'Location','NorthEast');
            legend boxoff
            grid on
        else
            text(.5,1.05,'Not computed','Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Cluster space
    axes('Units','pixels','Position',[840,130,255,170])
        fdata = mtint.clu_quality(sdata.tet(srow)).fetdata;
        nfets = mtint.clu_quality(sdata.tet(srow)).nfeatures;
        clus_cut = mtint.tetrode(sdata.tet(srow)).cut;
        clusters_now = unique(clus_cut);
        clus_count = numel(clusters_now);
        tspikes = mtint.clu_quality(sdata.tet(srow)).nspikes;
    
        % find the 2 channels with the maximum waveforms
%         if config.waveform_analyses
            [~,dindx] = sort(maxwavs,2,'descend'); % find the order of these values
%         else
%             dindx = [1 2];         
%         end
        mch1 = dindx(1); % the channel with the maximum amplitude
        mch2 = dindx(2); % the channel with the second maximum amplitude                    
        start_cols = 1:nfets:nfets*4; % the starting column for each channel    

        % get the feature data for the two channels
        d1 = fdata(:,start_cols(mch1)); % get this feature data for this channel
        d2 = fdata(:,start_cols(mch2)); % get this feature data for this channel

        % plot the noise in grey and then the cell in red
        colplot = [0.5 0.5 0.5 0.5];
        plot(d1,d2,'.','MarkerSize',3,'color',colplot); hold on; % plot the clusters
        cindx = find(clus_cut == sdata.clu(srow));
        plot(d1(cindx),d2(cindx),'.','MarkerSize',3,'color','r'); % re-plot the current cluster in red to make sure it is on top and stands out
        text(.5,1.05,sprintf('Ch%d vs Ch%d: %d clusters, %d spikes (%.1f%% of total %d)',mch1,mch2,clus_count,numel(pspx),numel(pspx)/tspikes*100,tspikes),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        xlabel(sprintf('PC1 Ch%d',mch1))
        ylabel(sprintf('PC1 Ch%d',mch2))
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
        rr = refline(1,0);
        set(rr,'Color',[.5 .5 .5],'LineStyle','-');
        uistack(rr,'bottom')
        grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Head direction (linear)
    axhd2 = axes('Units','pixels','Position',[1140,360,200,170]);   
        ai = rad2deg(linspace(0,2*pi,pdata.config.hd_bins)'); % angles for binning   
        hd_ratemap = sdata.hd_ratemap{srow};
        hd_dwellmap = pdata.(part_now).hd_dwellmap;

        % plot the data as a linear (side on) graph
        yyaxis left
        a1 = area(ai,hd_dwellmap,'FaceColor','k','FaceAlpha',.25);
        ylabel('Dwell time (s)')   
        set(gca,'yticklabel',num2str(get(gca,'ytick')','%.f'))                    
        yyaxis right
        a2 = area(ai,hd_ratemap,'FaceColor','b','FaceAlpha',.25);
        ylabel('Firing Rate (Hz)')        
        set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))            
   
        ax = gca;
        xlabel(sprintf('Head Direction (%s)',char(176)))
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'b';
        ax.XLim = [min(ai) max(ai)];
        ax.XTick = [0:90:360];
        text(.5,1.05,sprintf('r: %.2f, Max: %.2f, Mean: %.2f, Std: %.2f',sdata.hd_rayleigh(srow),sdata.hd_max(srow),sdata.hd_mean(srow),sdata.hd_stdev(srow)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Speed firing relationship  
    axvs = axes('Units','pixels','Position',[1140,130,200,170]);
        xi = 0:1:50;
        sf = sdata.speed_frate_curve{srow};
        sslope = sdata.speed_slope(srow);
        sscore = sdata.speed_score(srow);
        sitcpt = sdata.speed_y_intercept(srow);

        scatter(xi,sf,25,'ko','filled','MarkerFaceAlpha',0.5);
        ax = gca;
        ax.YLim(1) = 0;
        ax.XLim = [0 50];
        rr = refline(sslope,sitcpt);
        set(rr,'Color','r');
        
        xlabel('Speed (cm/s)')
        ylabel('Firing Rate (Hz)')
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k');
        box on; grid on;
        text(.5,1.05,sprintf('Score: %.2f, Slope: %.2f, Intercept: %.2f',sscore,sslope,sitcpt),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Overdispersion curve     
    axov = axes('Units','pixels','Position',[1390,590,200,170]);
        overd = sdata.over_dispersion(srow);  
        overr = sdata.over_dispersion_r(srow);          
        zdist = sdata.over_dispersion_z{srow};

        edg = -10:0.1:10;
        p1 = plot(edg,1,'k'); hold on;   
        if ~all(isnan(zdist(:)))
            delete(p1);
            fz = ksdensity(zdist,edg,'kernel','epanechnikov','Function','pdf');
            p1 = plot(edg,fz,'k');
        end

        % plot the expected gaussian
        fg = normpdf(edg,0,1);
        fg = fg ./ nanmax(fg(:)) .* nanmax(fz(:));
        p2 = plot(edg,fg,'r');

        ax = gca;
        ax.XLim = [-10 10];  
        xlabel('Standardized rate')
        ylabel('Probability')
        box on
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);         
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
        text(.5,1.05,sprintf('Mean: %.2f, OverD: %.2f, OverR: %.2f, N passes: %d',nanmean(zdist(:)),overd,overr,numel(zdist)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
        legend([p1 p2],{'Observed','Expected'},'Location','NorthEast');
        legend boxoff
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Cell ID text   
    axtxt = axes('Units','pixels','Position',[1390,120,200,170]);
        % get the cell's firing statistics
        frate = sdata.frate(srow);
        wwide = sdata.waveform_params(srow,2);
        skagg = sdata.spatial_info_bsec(srow);
        gscore = sdata.grid_score(srow);
        rvect = sdata.hd_rayleigh(srow);

        % place text in axis
        astr = sprintf('Mean firing rate: %.1f Hz\nWaveform width: %.1f ms\nSkaggs spatial information: %.1f b/s\nGrid score: %.1f\nRayleigh vector: %.1f\nAuto classified as: %s',frate,wwide,skagg,gscore,rvect,sdata.cell_type{srow});
        text(0.1,0.5,astr,'Units','normalized','FontSize',8,'FontWeight','normal')
        axis off    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Head direction polar plot 
    % WARNING: for some reason this plot must be done last, something in the mmpolar code messes up the axis calling
    axhd = axes('Units','pixels','Position',[1140,580,200,180]);
        set(fig_clust,'CurrentAxes',axhd);
        ai = deg2rad(ai);
        hd_dwellmap = hd_dwellmap ./ max(hd_dwellmap) .* max(hd_ratemap);       
        mmp = mmpolar(ai,hd_dwellmap,'k',ai,hd_ratemap,'b','FontSize',fsiz,'Grid','on','RGridVisible','off','RTickVisible','off','TTickDelta',20,'RTickLabelVisible','on','TTickLabelVisible','on');
        set(mmp(1),'LineWidth',0.5);
        set(mmp(2),'LineWidth',0.5);
        p1 = patch(get(mmp(1),'XData'),get(mmp(1),'YData'),'k','FaceAlpha',0.2);
        p2 = patch(get(mmp(2),'XData'),get(mmp(2),'YData'),'b','FaceAlpha',0.5);
        
        ax1 = get(gca,'Position');
        ll = legend([p1 p2],{'Session','Cell'},'Units','normalized','Position',[0.77 0.61 0 0],'FontSize',6,'NumColumns',2); legend boxoff;
        set(gca,'Position',ax1);
        legend boxoff
        axis square
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    figfile = [pwd '\klustest\' pdata.combined_name '\figures\'];
    [~,~,~] = mkdir(figfile);

    % saveas is faster than print
    saveas(gcf,[figfile num2str(sdata.rat(srow)) '_' num2str(sdata.date(srow)) '_' num2str(sdata.tet(srow)) '_' num2str(sdata.clu(srow)) '_' part_now '.png'],'png')
    close(fig_clust);                                          

















































