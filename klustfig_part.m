function klustfig_part(pdata,sdata,pp,wav_now,quals,fets,fast_figs,sshuff,figvis)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  short description
% long description
%
% USAGE:
%       [out] = template(in,in2)
%
% INPUT:
%       in - input 1
%       in2 - input 2
%
% OUTPUT:
%       p - output
%
% EXAMPLES:
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 13/02/23 Initial release
% version 1.1.0, Release 13/02/23 Spike plot, firing rate map, autocorrelation working
% version 1.1.1, Release 13/02/23 Using new 'turbo' colormap for grid autocorrelation
% version 1.1.2, Release 13/02/23 Using 'viridis' colormap for firing rate map
% version 1.2.0, Release 14/02/23 Waveform plots working
% version 1.3.0, Release 15/02/23 Head direction, AHV, speed plots working, simplified HD polar plot
% version 1.4.0, Release 15/02/23 Spike autocorrelograms working, added theta fit plot below 500ms autocorr, ISI plot working
% version 1.4.1, Release 15/02/23 Spikes within refractory period are colored red
% version 1.5.0, Release 15/02/23 Spike histogram along bottom added
% version 1.6.0, Release 16/02/23 Feature space & cluster quality plots working, CDFs added to cluster quality
% version 1.6.1, Release 16/02/23 Annotation & uci added to top of plot
% version 2.0.0, Release 16/02/23 Added faster figure saving method using getframe
% version 3.0.0, Release 18/02/23 Added limits to the number of waveforms and features plotted to prevent crashing
% version 3.0.0, Release 18/02/23 Added text to show the N spikes and features displayed
% version 3.0.0, Release 18/02/23 Using plot for spikes instead of scatter, faster and looks better I think
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PREPARE FIGURE
    fig_clust = figure('visible',figvis,'Units','pixels','Position',[10, 10, 1600, 900]); % open/create figure
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    fsiz = 8; % the fontsize for different texts
    flnw = 0.5; % the line width for different plots
% set(gcf,'visible','on')
    % add an annotation to the figure with some important info
    part_now = pdata.part_config.part_names{pp};
    tet_now = sdata.tetrode(1);
    clu_now = sdata.cluster(1);    
    ann_str = sprintf('Cell: %s, Part: %s, Analysed: %s',sdata.uci{1},part_now,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
    annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',8,'LineStyle','none','interpreter','none');  
    
    % settings
    max_plot_waves = 250; % set to 0 to show all waveforms
    max_plot_fets = 1e5; % set to 0 to show all feature points in feature space
    max_plot_spikes = 1e5; % set to 0 to show all spikes in spike & position plot
    tets = pdata.tetrodes;
    mapset = pdata.mapset;
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Positions and spikes
    axps = axes('Units','pixels','Position',[30,645,400,200]);
        % position and spike data (already cut to this part)
        part_duration = pdata.part_config.part_duration(pp);
        ppox = pdata.(part_now).pox;
        ppoy = pdata.(part_now).poy;
        ppot = pdata.(part_now).pot; 
        pspx = ppox(sdata.spt_pot_index{1});
        pspy = ppoy(sdata.spt_pot_index{1});
        
        % plot position data, excluding pieces not included in this part
        % by inserting NaNs between intervals we can plot this as one line, which saves on memory
        dindax = abs([0; diff(ppot)])>0.1;
        pos_plot = [ppox ppoy];
        pos_plot(dindax,:) = NaN;
        plot(pos_plot(:,1),pos_plot(:,2),'Color',[.5 .5 .5 .5]); hold on;

        % plot spikes after position data so they are all on top
        if max_plot_spikes>0 % if we want to limit the number of waveforms shown
            max_plot_spikes_now = min([numel(pspx) max_plot_spikes]); % if there are less than max_plot_spikes, plot them all
            rindx = randperm(numel(pspx),max_plot_spikes_now);
            spx_plot = pspx(rindx,:);
            spy_plot = pspy(rindx,:); 
        else
            spx_plot = pspx;
            spy_plot = pspy;             
        end        
        plot(spx_plot,spy_plot,'Color',[1 0 0 0.5],'Marker','.','LineStyle','none') 
        
        % additional settings
        if mapset.fix_aspect % if we want the plot to fill the figure (i.e. plot in landscape if the data are landscape)
            a_ratio = range(ppox) / range(ppoy);
            if a_ratio < 1 % y spans more distance than x
                view(90,90); % rotate the plot
            else
                view(0,90);
            end
        end
        daspect([1 1 1])
        axis xy off tight
        text(0,1.1,sprintf('%d spikes (%.2f Hz), %d seconds (%.1f mins)',numel(pspx),numel(pspx)/part_duration,round(part_duration),part_duration/60),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
        
%% >>>>>>>>>> Firing rate map       
    axrt = axes('Units','pixels','Position',[axps.Position(1) axps.Position(2)-axps.Position(4)-55 axps.Position(3) axps.Position(4)]);
        ratemap = sdata.ratemap{1};
        spati = sdata.spatial_info(1);
        spars = sdata.spatial_info(3);
        coher = sdata.spatial_info(4);

        im = imagesc(ratemap,'alphadata',~isnan(ratemap));
        daspect([1 1 1])
        clim(double([0 max([0.1 max(ratemap(:),[],'omitnan')])]))       
        colormap(axrt,turbo);
        axis xy off
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);   
        if mapset.fix_aspect % if we want the plot to fill the figure (i.e. plot in landscape if the data are landscape)
            if a_ratio < 1 % y spans more distance than x
                view(90,90); % rotate the plot
            else
                view(0,90);
            end
        end

        if sshuff % spatial shuffles were performed
            spatiz = sdata.spatial_info_z(1);
            spatip = sdata.spatial_info_p(1);
            text(0,1.05,sprintf('SI: %.2f (z %.2f,p %.3f), Sp: %.2f, Cohe: %.2f',spati,spatiz,spatip,spars,coher),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
        else
            text(0,1.05,sprintf('SI: %.2f, Sp: %.2f, Cohe: %.2f',spati,spars,coher),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
        end

        axp = get(gca,'Position');
        cc = colorbar;
        set(gca,'Position',axp);
        set(cc,'Position',get(cc,'Position')+[0 0 -0.004 0])
        title(cc,'Hz','FontSize',fsiz)        
        set(cc,'yticklabel',num2str(get(cc,'ytick')','%.1f')) % change Ytick labels to have the same number of decimal places
        
%% >>>>>>>>>> Grid autocorrelogram              
    axgs = axes('Units','pixels','Position',[axrt.Position(1) axrt.Position(2)-axrt.Position(4)-55 axrt.Position(3) axrt.Position(4)]);
        amap = sdata.gridmap{1};  
        gscore = sdata.grid_info(1);
        gspacing = sdata.grid_info(2);
        gori = sdata.grid_info(3);        
        
        imc = imagesc(amap,'alphadata',~isnan(amap));
        daspect([1 1 1])
        clim(double([-0.2 1]))
        colormap(axgs,hot);        
        axis xy off
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
        if mapset.fix_aspect % if we want the plot to fill the figure (i.e. plot in landscape if the data are landscape)
            if a_ratio < 1 % y spans more distance than x
                view(90,90); % rotate the plot
            else
                view(0,90);
            end
        end      

        if sshuff % spatial shuffles were performed
            gridz = sdata.spatial_info_z(2);
            gridp = sdata.spatial_info_p(2);
            text(0,1.05,sprintf('G: %.2f (z %.2f,p %.3f), W: %.2f, O: %.2f',gscore,gridz,gridp,gspacing,gori),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')        
        else
            text(0,1.05,sprintf('G: %.2f, W: %.2f, O: %.2f',gscore,gspacing,gori),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')        
        end

        axp = get(gca,'Position');
        cc = colorbar;
        set(gca,'Position',axp);
        set(cc,'Position',get(cc,'Position')+[0 0 -0.004 0])
        title(cc,'Hz','FontSize',fsiz)        
        set(cc,'yticklabel',num2str(get(cc,'ytick')','%.1f')) % change Ytick labels to have the same number of decimal places
                
%% >>>>>>>>>> Waveforms            
    axwav = axes('Units','pixels','Position',[axps.Position(1)+axps.Position(3)+80 axps.Position(2) 200 axps.Position(4)]);
        mxs = NaN(length(wav_now),1);
        wav_means = NaN(length(wav_now),size(wav_now{1},2));
        wav_stds = NaN(length(wav_now),size(wav_now{1},2));  
        ax_lim = [-0.1 0.1];
        for ww = 1:length(wav_now) % for every channel
            wnow = wav_now{ww}; % should be [spikes x samples]
            if isempty(wnow)
                continue
            end
            wav_means(ww,:) = mean(wnow,1,'omitnan');
            wav_stds(ww,:) = std(double(wnow),[],1,'omitnan');
            mxs(ww,1) = max(wav_means(ww,:),[],'omitnan');
%             ax_lim(1) = min([ax_lim(1) min(wnow(:))]);
%             ax_lim(2) = max([ax_lim(2) max(wnow(:))]);
            ax_lim(1) = min([ax_lim(1) min(wav_means(ww,:)-(wav_stds(ww,:).*3))]);
            ax_lim(2) = max([ax_lim(2) max(wav_means(ww,:)+(wav_stds(ww,:).*3))]);
        end
        [~,widx] = sort(mxs,'descend','MissingPlacement','last'); % sort from largest > smallest waveform      

        % plot the waveform(s) with the highest amplitude
        wnow = squeeze(wav_now{widx(1)});
        if max_plot_waves>0 % if we want to limit the number of waveforms shown
            max_plot_waves_now = min([size(wnow,1) max_plot_waves]); % if there are less than max_plot_waves waveforms, plot them all
            rindx = randperm(size(wnow,1),max_plot_waves_now);
            wnow = wnow(rindx,:);
        end
        
        wavtime = pdata.wavtime; 
        if isempty(wnow) % if there are no spikes
            wnow = NaN(size(wavtime'));
        end

        plot(wavtime,wnow','Color','k'); hold on;
        [hl,hp] = boundedline(wavtime,wav_means(widx(1),:),wav_stds(widx(1),:),'-k');
        set(hl,'Color','r','LineStyle','-','LineWidth',1); 
        set(hp,'FaceColor','b','FaceAlpha',.5);     
        ax = gca;
        
        ax.XLim = mapset.wave_window;
        ax.YLim = ax_lim;
        ax.YDir = 'normal';
        box on
        grid on
        xlabel('Time (ms)')
        ylabel(sprintf('Amplitude (%cV)',char(956)))
        text(1,1.05,sprintf('Ch%d, Peak: %.1f%cV, Width: %.2fms',widx(1),mxs(widx(1)),char(956),sdata.wave_width),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','right')
        if max_plot_waves>0 % if we want to limit the number of waveforms shown
            text(0.99,1,sprintf('N = %d',max_plot_waves_now),'Units','normalized','FontSize',8,'HorizontalAlignment','right','VerticalAl','top')           
        end

        % plot all 4 waveforms in smaller plots
        siz = 50;
        xbuff = 20;
        ybuff = siz;
        for ww = 1:4
            axes('Units','pixels','Position',[axwav.Position(1)+axwav.Position(3)+xbuff axwav.Position(2)+axwav.Position(4)-(ww*ybuff) siz siz])
                wnow = squeeze(wav_now{ww});
                
                if max_plot_waves>0 % if we want to limit the number of waveforms shown
                    max_plot_waves_now = min([size(wnow,1) max_plot_waves]); % if there are less than max_plot_waves waveforms, plot them all
                    rindx = randperm(size(wnow,1),max_plot_waves_now);
                    wnow = wnow(rindx,:);
                end   
                
                if isempty(wnow) % if there are no spikes
                    wnow = NaN(size(wavtime'));
                end
                plot(wavtime,wnow','Color','k'); hold on;
                [hl,hp] = boundedline(wavtime,wav_means(ww,:),wav_stds(ww,:),'-k');
                set(hl,'Color','r','LineStyle','-','LineWidth',1); 
                set(hp,'FaceColor','b','FaceAlpha',.5);  
                ax = gca;
                ax.XLim = mapset.wave_window;
                ax.YLim = ax_lim;
                ax.XTick = [];
                ax.YTick = [];
                ax.YDir = 'normal';                
                box on   
                grid on
                set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);   
                text(-0.2,0.5,sprintf('Ch%d',ww),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center','rotation',90)                
        end          
        
%% >>>>>>>>>> Isolation quality
        axiso = axes('Units','pixels','Position',[axwav.Position(1)+20 axwav.Position(2)-axwav.Position(4)-55 axwav.Position(3)+30 axwav.Position(4)]);
            dists = double(quals{tet_now}{clu_now});
            if ~isempty(dists)
                h1 = dists(2,:)./sum(dists(2,:));
                h2 = dists(3,:)./sum(dists(3,:));
                b1 = bar(dists(1,:),h1,1,'b','FaceAlpha',0.5,'EdgeColor','none'); hold on;
                b2 = bar(dists(1,:),h2,1,'k','FaceAlpha',0.5,'EdgeColor','none');         

                % plot inset CDF
                scaler = 0.4;
                ax_ins = axes('Units','pixels','Position',[axiso.Position(1)+(axiso.Position(3)-axiso.Position(3)*scaler) axiso.Position(2)+(axiso.Position(4)-axiso.Position(4)*scaler) axiso.Position(3)*scaler axiso.Position(4)*scaler]);
                    plot(dists(1,:),cumsum(h1),'Color','b'); hold on;
                    plot(dists(1,:),cumsum(h2),'Color','k');            
                    ax_ins.XScale = 'log';
                    ax_ins.XLim = axiso.XLim;
                    ax_ins.XTick = [];
                    ax_ins.YTick = [];
            end
            set(gcf,'currentaxes',axiso); % set current axis in figure
            axiso.XScale = 'log';
            xlabel('Mahalanobis distance')
            ylabel('Prop. distances')   
            set(axiso,'yticklabel',num2str(get(gca,'ytick')','%.2f'))                                
            axiso.XTick = [0 1 10 100 1000 10000 100000];
            text(0,1.05,sprintf('Iso-d: %.2f, Lratio: %.3f',sdata.isod(1),sdata.isod(2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')            
            
                                
%% >>>>>>>>>> Feature space
        axfet = axes('Units','pixels','Position',[axiso.Position(1) axiso.Position(2)-axiso.Position(4)-55 axiso.Position(3) axiso.Position(4)]);
            clu = pdata.clusters{tet_now};
            fet = fets{tet_now};
            maxwavs = max(sdata.wave_mean{1},[],2);
            [~,dindx] = sort(maxwavs,1,'descend'); % find the order of these values
            nch = sum(~cell2mat(cellfun(@isempty,wav_now,'UniformOutput',false)));

            nfets = (size(fet,2))/nch;            
            start_cols = 1:nfets:nfets*nch; % the starting column for each channel    
            fet_to_plot = 4;
            % get the feature data for the two max channels
            d1 = fet(:,start_cols(dindx(1))+(fet_to_plot-1)); % get this feature data for this channel
            d2 = fet(:,start_cols(dindx(2))+(fet_to_plot-1)); % get this feature data for this channel

            % downsample feature space
            if max_plot_fets>0 % if we want to limit the number of feature space points shown
                max_plot_fets_now = min([size(d1,1) max_plot_fets]); % if there are less than max_plot_fets points, plot them all
                rindx = randi(size(d1,1),[max_plot_fets_now,1]);
                d1 = d1(rindx,:);
                d2 = d2(rindx,:);
                clu = clu(rindx);
            end            
            
            % plot the clusters    
            cols = winter(double(max(clu)));
            for cc = 1:max(clu)
                if cc==clu_now % skip the current cluster because we want to do that last & in black
                    continue
                end
                plot(d1(clu==cc),d2(clu==cc),'.','MarkerSize',4,'color',cols(cc,:)); hold on; % plot the clusters
            end
            plot(d1(clu==clu_now),d2(clu==clu_now),'.','MarkerSize',6,'color','k'); % re-plot the current cluster in red to make sure it is on top and stands out

            axis tight
            xlabel(sprintf('Fet1 Ch%d',dindx(1)))
            ylabel(sprintf('Fet1 Ch%d',dindx(2)))              
            grid on            
            if max_plot_fets>0 % if we want to limit the number of waveforms shown
                text(0.99,1,sprintf('N = %d',max_plot_fets_now),'Units','normalized','FontSize',8,'HorizontalAlignment','right','VerticalAl','top')           
            end
        
%% >>>>>>>>>> Head direction
        axhd = axes('Units','pixels','Position',[axwav.Position(1)+axwav.Position(3)+120 axwav.Position(2)+45 axwav.Position(3) axwav.Position(4)-45]);
            ai = rad2deg(linspace(0,2*pi,pdata.mapset.hd_bins)'); % angles for binning   
            hd_ratemap = sdata.hd_ratemap{1};
            hd_dwellmap = pdata.(part_now).hd_dwellmap;

            % plot the data as a linear (side on) graph
            a2 = area(ai,hd_ratemap,'FaceColor','k','FaceAlpha',.75,'EdgeColor','none'); hold on;
            ylabel('Firing Rate (Hz)')        

            ax = gca;
            ax.XLim = [min(ai) max(ai)];
            ax.XTick = [0:90:360]; 
            ax.XTickLabel = {};

            if sshuff % spatial shuffles were performed
                rz = sdata.spatial_info_z(3);
                rp = sdata.spatial_info_p(3);
                text(0,1.05,sprintf('r: %.2f (z %.2f,p %.3f), PFD%c: %.2f',sdata.hd_info(1),rz,rp,176,sdata.hd_info(2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
            else
                text(0,1.05,sprintf('r: %.2f, PFD%c: %.2f',sdata.hd_info(1),176,sdata.hd_info(2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
            end

            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);       
            grid on;
            ax.YLim = [0 max([0.1 max(hd_ratemap(:),[],'omitnan')]) ].*1.1;
            line([sdata.hd_info(2) sdata.hd_info(2)],ax.YLim,'Color','r')
            set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))   
       
        axhd2 = axes('Units','pixels','Position',[axhd.Position(1) axhd.Position(2)-45 axhd.Position(3) 40]);
            a1 = area(ai,hd_dwellmap,'FaceColor','b','FaceAlpha',.75,'EdgeColor','none'); hold on;
            ylabel('Time (s)')   
            set(gca,'yticklabel',num2str(get(gca,'ytick')','%.f'))            
        
            ax = gca;
            ax.YAxisLocation = 'right';
            xlabel(sprintf('Head Direction (%s)',char(176)))
            ax.XLim = [min(ai) max(ai)];
            ax.XTick = [0:90:360];        
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);           
            line([sdata.hd_info(2) sdata.hd_info(2)],ax.YLim,'Color','r')
            grid on;
        
%% >>>>>>>>>> Head direction polar plot
        axhd2 = polaraxes('Units','pixels','Position',[axhd.Position(1) axhd.Position(2)-axhd.Position(4)-155 axhd.Position(3) axhd.Position(4)+25]);
            theta = linspace(0,2*pi,pdata.mapset.hd_bins)'; 
            p1 = polarplot(theta,hd_ratemap,'k'); hold on;
            dwellmap_to_plot = (hd_dwellmap/max(hd_dwellmap(:)))*max(hd_ratemap)*0.5;
            p1 = polarplot(theta,dwellmap_to_plot,'b'); hold on;            
            axhd2.ThetaZeroLocation = 'right';

            fname = 'C:\Users\F004KS7\Downloads\image.png';
            exportgraphics(gcf,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  

%% >>>>>>>>>> HD half plots
        axav = axes('Units','pixels','Position',[axhd2.Position(1) axhd2.Position(2)-axhd2.Position(4)-70 axhd2.Position(3) axhd2.Position(4)]);
            ai = rad2deg(linspace(0,2*pi,pdata.mapset.hd_bins)'); % angles for binning   
            hd_ratemap1 = sdata.hd_ratemap_half{1};
            hd_ratemap2 = sdata.hd_ratemap_half{2};

            % plot the data as a linear (side on) graph
            a1 = area(ai,hd_ratemap1,'FaceColor','k','FaceAlpha',.75,'EdgeColor','none'); hold on;
            a2 = area(ai,hd_ratemap2,'FaceColor','b','FaceAlpha',.75,'EdgeColor','none'); hold on;            
            ylabel('Firing Rate (Hz)')        

            ax = gca;
            ax.XLim = [min(ai) max(ai)];
            ax.XTick = [0:90:360]; 
            ax.XTickLabel = {};
            text(0,1.05,sprintf('r: %.2f, PFD%c: %.2f | r: %.2f, PFD%c: %.2f',sdata.hd_info_half(1),176,sdata.hd_info_half(2),sdata.hd_info_half(5),176,sdata.hd_info_half(6)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
            ax.YLim = [0 max([0.1 max(hd_ratemap1(:),[],'omitnan') max(hd_ratemap2(:),[],'omitnan')]) ].*1.1;
            line([sdata.hd_info_half(2) sdata.hd_info_half(2)],ax.YLim,'Color','k')
            line([sdata.hd_info_half(6) sdata.hd_info_half(6)],ax.YLim,'Color','b')
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);       
            grid on;
            [~,leg] = legendflex([a1 a2],{'Half 1','Half 2'},'anchor',{'nw','nw'},'ncol',1,'box','off','buffer',[0,0],'fontsize',fsiz);    
            leg(3).Children.FaceAlpha = 0.75;
            leg(4).Children.FaceAlpha = 0.75;
            set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))   

%% >>>>>>>>>> ISIs
        axisi = axes('Units','pixels','Position',[axhd.Position(1)+axhd.Position(3)+70 axhd.Position(2)-45 axhd.Position(3) axhd.Position(4)+45]);
            xi = [-pdata.isi_xvalues pdata.isi_xvalues];
            idist = [sdata.isi{1}; sdata.isi{1}];
            bar(xi(:),idist(:),1,'k','EdgeColor','k'); hold on;
            bar(xi(abs(xi)<2),idist(abs(xi)<2),1,'r','EdgeColor','r');

            ax = gca;
            ax.XLim = [-50 50];
            ax.XTick = -50:25:50;                        
            ylabel('Spikes')
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); 
            xlabel('Interspike interval (ms)') % label x axis
            text(.5,1.05,sprintf('fwhmx: %.1fms, burst index: %.2f',sdata.isi_info(1),sdata.burst_index(1)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
            grid on;
            %line([-2 -2],ax.YLim,'Color','b')
            %line([2 2],ax.YLim,'Color','b')
        
%% >>>>>>>>>> Spike autocorrelation - 25ms refractory period             
        axref = axes('Units','pixels','Position',[axisi.Position(1) axisi.Position(2)-265 axisi.Position(3) axisi.Position(4)-10]);
            autoc = sdata.autocorr_25{1};
            tlag = pdata.autocorr_25_xvalues;

            bar(tlag,autoc,1,'k','EdgeColor','k'); hold on;
            bar(tlag(abs(tlag)<2),autoc(abs(tlag)<2),1,'r','EdgeColor','r'); hold on;
            
            ax = gca;            
            ax.XLim = [-20 20];
            ax.XTick = -20:10:20;            
            xlabel('Time lag (ms)') % label x axis
            ylabel('Probability') % label y axis
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
            text(0,1.05,sprintf('RPV: %d (%.2f%%)',sdata.autocorr_25_info(1),sdata.autocorr_25_info(2)*100),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
            set(ax,'yticklabel',num2str(get(gca,'ytick')','%.3f'))        
            grid on;
            %line([-2 -2],ax.YLim,'Color','b')
            %line([2 2],ax.YLim,'Color','b')
            
%% >>>>>>>>>> Spike autocorrelation - 500ms theta modulation           
        axref = axes('Units','pixels','Position',[axref.Position(1) axref.Position(2)-210 axref.Position(3) axref.Position(4)-40]);
            autoc = sdata.autocorr_500{1}(:,1);
            tlag = pdata.autocorr_500_xvalues;
            autof = sdata.autocorr_500{1}(:,2);

            bar(tlag,autoc,1,'k'); hold on;

            ax = gca;
            ax.XLim = [-500 500];
            ax.XTick = -500:125:500;
            ax.XTickLabel = {};
            ylabel('Probability') % label y axis
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
            text(.5,1.05,sprintf('Theta index: %.1f, frequency: %.2fHz',sdata.autocorr_500_info(1),sdata.autocorr_500_info(2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
            set(ax,'yticklabel',num2str(get(gca,'ytick')','%.3f'))        
            grid on;
            
        axref2 = axes('Units','pixels','Position',[axref.Position(1) axref.Position(2)-45 axref.Position(3) 40]);
            a2 = area(tlag,autof,'FaceColor','b','FaceAlpha',.75,'EdgeColor','none'); hold on;
            ax = gca;
            ax.XLim = [-500 500];
            ax.XTick = -500:125:500;
            ax.XTickLabel = {'-500','','-250','','0','','250','','500'};
            ax.XTickLabelRotation = 0;
            ax.YAxisLocation = 'right';            
            xlabel('Time lag (ms)') % label x axis
            set(ax,'yticklabel',num2str(get(gca,'ytick')','%.2f'))        
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
            line([0 0],ax.YLim,'Color',[.5 .5 .5])
            grid on
            
%% >>>>>>>>>> Theta phase preference
        axthp = axes('Units','pixels','Position',[axisi.Position(1)+axisi.Position(3)+60 axisi.Position(2) axisi.Position(3) axisi.Position(4)]);
            phase_dist = sdata.theta_phase{1};
            ai = reshape(deg2rad(-180:5:540),[],1); % doing this means we have bins symmetrical around zero
            xi = movmean(ai,2,'EndPoints','discard');

            bar(xi,phase_dist,1,'k'); hold on;
            ax = gca;
            y1 = ax.YLim;
            si = cos(xi) ./ max(sin(xi)) .* (y1(2)/2) + (y1(2)/2);
            plot(xi,si,'Color','r','LineWidth',1);
            ax.YLim = y1;

            ylabel('Spikes')
            ax.XTick = [-pi 0 pi 2*pi 3*pi];
            ax.XTickLabel = {'-\pi','0','\pi','2\pi','3\pi'};
            xlabel('Theta phase (rad)') % label x axis
            text(.5,1.05,sprintf('R: %.1f, PFA: %.2f',sdata.theta_info(1),sdata.theta_info(2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')    
            
%% >>>>>>>>>> Speed
        axvs = axes('Units','pixels','Position',[axthp.Position(1) axthp.Position(2)-265 axthp.Position(3) axthp.Position(4)-10]);
            xi = 0:1:50;
            sf = sdata.speed_slope{1};
            sscore = sdata.speed_info(1);
            sslope = sdata.speed_info(2);
            sitcpt = sdata.speed_info(3);

            scatter(xi,sf,25,'ko','filled','MarkerFaceAlpha',0.5);
            ax = gca;
            ax.YLim(1) = 0;
            ax.XLim = [0 50];
            rr = refline(sslope,sitcpt);
            set(rr,'Color','r');

            xlabel('Speed (cm/s)')
            ylabel('Firing Rate (Hz)')
            set(ax,'yticklabel',num2str(get(gca,'ytick')','%.1f'))                    
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k');
            box on; grid on;
            text(0,1.05,sprintf('Score: %.2f, Slope: %.2f, Intercept: %.2f',sscore,sslope,sitcpt),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')        
                    
%% >>>>>>>>>> AHV
        axav = axes('Units','pixels','Position',[axvs.Position(1) axvs.Position(2)-axvs.Position(4)-70 axvs.Position(3) axvs.Position(4)]);
            xi = pdata.ahv_xvalues;
            ahv = sdata.ahv_curve{1};
            ahv(ahv==0) = NaN; % don't show zero values
            scatter(xi,ahv,25,'ko','filled','MarkerFaceAlpha',0.5);

            ax = gca;
            ax.YLim(1) = 0;
            xspan = 250; % deg/s
            ax.YLim(2) = max([0.001 max(ahv(abs(xi)<xspan))*1.1]);
            ax.XLim = [-xspan xspan];     
            xlabel(sprintf('AHV (%c/s)',176))
            ylabel('Firing Rate (Hz)')
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k');
            line([0 0],ax.YLim,'Color',[.5 .5 .5]);
            box on; grid on;

%% >>>>>>>>>> Spike rasterplot along bottom of figure
    axsr = axes('Units','normalized','Position',[0.025,0.045,0.945,0.045]);
        % get the spike time data and bin it
        bsize = 5; 
        pos = pdata.pos;
        edg = min(pos.pot) : bsize : max(pos.pot); % vector of time points at which we should calculate spike frequency
        spt = pdata.spike_times{tet_now}(pdata.clusters{tet_now}==clu_now); % all spike times for this cluster (including outside the current part)
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
        tvals = pdata.part_config.part_times{pp};

        for tt = 1:length(tvals(:,1))
            a1 = line([knnsearch(edg(:),tvals(tt,1)) knnsearch(edg(:),tvals(tt,2))],[0.5 0.5],'Color','r','LineWidth',4);
        end
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % set the axes to be on top, set their thickness to flnw, set all fonts to size fsiz and make sure the axis and text are black 
        set(ax,'xticklabel',num2str(get(gca,'xtick')'.*bsize,'%.f'))
        %[~,leg] = legendflex([a1],{'= current part'},'anchor',{'sw','sw'},'ncol',1,'box','off','buffer',[-20,-35],'fontsize',fsiz);    

%% >>>>>>>>>> Save the overall figure
        % Save the figure  
        [~,~,~] = mkdir([pwd '\' pdata.outname '\part_figures']); % create a folder to hold outputs 
        fname = [pwd '\' pdata.outname '\part_figures\' sdata.uci{1} '_' part_now '.png'];
        warning('off','MATLAB:Axes:NegativeLimitsInLogAxis');
        if fast_figs
            frame = getframe(fig_clust); % fig is the figure handle to save
            [raster, raster_map] = frame2im(frame); % raster is the rasterized image, raster_map is the colormap
            if isempty(raster_map)
                imwrite(raster, fname);
            else
                imwrite(raster, raster_map, fname); % fig_file is the path to the image
            end            
        else
            exportgraphics(fig_clust,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',150);  
        end
        close(fig_clust);    




































