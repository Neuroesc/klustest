function figCLUS(mtint,pdata,sdatac,fig_vis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figCLUS  a basic figure function
% This function just stores the plotting for klustest. This figure gives all the information about a cluster
% in all parts (i.e. spikes, ratemap)
%
% USAGE:
%         figCLUS(mtint,pdata,sdata)
%
% INPUT:
%         mtint - mtint structure from klustest
%         pdata - pdata structure from klustest
%         sdata - sdatap structure from klustest
%
% See also: klustest figPART

% HISTORY:
% version 1.0.0, Release 12/04/19 Initial release, created to replace what was figPARTS
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with input variables
inps = {'fig_vis','save_fig'};
vals = {'''off''','0'};
for ff = 1:length(inps)
    if ~exist(inps{ff},'var')
        eval([inps{ff} '=' vals{ff} ';']);
    end
end

if ~any(sdatac.spikes) || all(isempty(sdatac.spikes)) || all(isnan(sdatac.spikes))
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
    ann_str = sprintf('Cell: %s, Rat: %s, Date: %d, Tetrode: %d, Cluster: %d, Analysed: %s, Auto cell classification: %s',sdatac.uci{1},sdatac.rat{1},sdatac.date(1),sdatac.tet(1),sdatac.clu(1),datestr(now,'yyyy-mm-dd-HH-MM-SS'),sdatac.cell_type{1});
    annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',8,'LineStyle','none','interpreter','none');  

    % collect data before parts loop
    part_names = sdatac.part(:);
    
    % subplot controls
    fh = length(part_names); % how many plots wide should it be
    fv = 2; % how many plots tall should it be
    fs = 0.03; % the spacing around the plots, on all sides
    fp = 0.02; % the spacing around the plots, on all sides, this takes more space than fspac though
    fm = 0.08; % the margins around the plots, at the edge of the figure  
    ax_mat = reshape((1:(fh*fv)),fh,fv)';
    
    for pp = 1:length(part_names) % for every part in which the cell was recorded
        part_now = part_names{pp};
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Position and spike plot (black lines & red dots)    
        subaxis(fv,fh,ax_mat(1,pp),'Spacing',fs,'Padding',fp,'Margin',fm,'Holdaxis');
            % position and spike data (already cut to this part)
            part_duration = sdatac.duration(pp,1);  
            ppox = pdata.(part_now).pox;
            ppoy = pdata.(part_now).poy;
            ppot = pdata.(part_now).pot; 
            pspx = pdata.pox(sdatac.spike_index{pp,1});
            pspy = pdata.poy(sdatac.spike_index{pp,1});

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
        
            if numel(pspx)<=pdata.minspikes
                continue
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Firing rate map
        axr = subaxis(fv,fh,ax_mat(2,pp),'Spacing',fs,'Padding',fp,'Margin',fm,'Holdaxis');
            ratemap = sdatac.ratemap{pp,1};
            spati = sdatac.spatial_info_bsec(pp,1);
            spars = sdatac.sparsity(pp,1);
            coher = sdatac.spatial_coherence(pp,1);

            im = imagesc(ratemap);
            set(im,'alphadata',~isnan(ratemap));
            daspect([1 1 1])
            caxis([0 nanmax(ratemap(:))])        
            axis xy off tight square
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);   

            text(.5,1.05,sprintf('SI: %.2f, Sp: %.2f, Cohe: %.2f',spati,spars,coher),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
            tt = text(-0.07,0.5,sprintf('Max: %.2f, Mean: %.2f, Median: %.2f',nanmax(ratemap(:)),nanmean(ratemap(:)),nanmedian(ratemap(:))),'FontSize',fsiz+1,'Units','normalized','HorizontalAlignment','center');
            set(tt,'rotation',90);

            colormap(axr,'jet');
            axp = get(gca,'Position');
            cc = colorbar;
            set(gca,'Position',axp);
            set(cc,'Position',get(cc,'Position')+[0 0 -0.004 0])
            title(cc,'Hz','FontSize',fsiz)        
            
    end % ends the parts loop
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Spike rasterplot running along bottom of figure
    axsr = axes('Units','pixels','Position',[30,35,1550,45]);
        % get the spike time data and bin it
        bsize = 5;
        edg = 0:bsize:pdata.duration; % vector of 1s time points at which we should calculate spike probability
        spt = mtint.tetrode(sdatac.tet(1,1)).ts( mtint.tetrode(sdatac.tet(1,1)).cut==sdatac.clu(1,1) );
        f = histcounts(spt,edg);

        % plot this as an image
        im = imagesc(f); hold on;
        colormap(axsr,flipud(gray(256)))
        ax = gca;
        ax.YTick = [];
        ax.TickLength = [0.001, 0.001];
        ylabel('Spikes')
        hold on
        axis xy on
        xlabel('Time (s)')
        
        cols = hsv(length(part_names));
        for pp = 1:length(part_names) % for every part in which the cell was recorded
            tvals = pdata.(part_names{pp}).times;
            for tt = 1:length(tvals(:,1))
                line([knnsearch(edg(:),tvals(tt,1)) knnsearch(edg(:),tvals(tt,2))],[0.5 0.5],'Color',cols(pp,:),'LineWidth',4);
            end
        end
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % set the axes to be on top, set their thickness to flnw, set all fonts to size fsiz and make sure the axis and text are black 
        set(ax,'xticklabel',num2str(get(gca,'xtick')'.*bsize,'%.f'))
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    figfile = [pwd '\klustest\' pdata.combined_name '\figures\'];
    [~,~,~] = mkdir(figfile);

    % saveas is faster than print
    saveas(gcf,[figfile sdata.uci{1} '_xsummary.png'],'png')
    close(fig_clust);                                          

















































