function klustfig_hills(pdata,cluma,uci,fast_figs,figvis,fname)
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
% version 1.0.0, Release 19/02/23 Initial release
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PREPARE FIGURE
    fig_clust = figure('visible',figvis,'Units','pixels','Position',[10, 10, 1500, 800]); % open/create figure
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    fsiz = 8; % the fontsize for different texts
    flnw = 0.5; % the line width for different plots

    % add an annotation to the figure with some important info
    ann_str = sprintf('Cell: %s, Analysed: %s',uci,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
    annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',8,'LineStyle','none','interpreter','none');  
    
    part_names = {'arena1','hills','arena2'};
    nparts = numel(unique(pdata.pos.session));   
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    xmin = 50;
    xbuff = 350;
    xvec = [xmin xmin+xbuff xmin+2*xbuff xmin+3*xbuff];
    ymax = 500;
    ybuff = 200;
    yvec = [ymax ymax-ybuff ymax-2*ybuff ymax-3*ybuff];
    pwidth = 300;
    pheight = 250;
    
    part_idx = [1 2 2 3];
    type_idx = {'planar','planar','surficial','planar'};
    for cc = 1:4 % for each column of plots
        part_now = part_names{part_idx(cc)};
        idx = find( ismember(cluma.uci,uci) & cluma.session==part_idx(cc) );
        pidx = pdata.pos.session==part_idx(cc);
        part_duration = cluma.session_duration_s(idx);
        ppox = pdata.pos.(['pox_' type_idx{cc}])(pidx);
        ppoy = pdata.pos.(['poy_' type_idx{cc}])(pidx);
        ratemap = cluma.(['ratemap_' type_idx{cc}]){idx};
        amap = cluma.([type_idx{cc} '_amap']){idx};         
        ppot = pdata.pos.pot(pidx);         
        pspx = pdata.pos.(['pox_' type_idx{cc}])(cluma.spike_index{idx});
        pspy = pdata.pos.(['poy_' type_idx{cc}])(cluma.spike_index{idx});  
        sinfo = cluma.([type_idx{cc} '_spatial_info_shuffles'])(idx,:);
        
%% >>>>>>>>>> Positions and spikes
    axps = axes('Units','pixels','Position',[xvec(cc),yvec(1),pwidth,pheight]);
        % plot position data, excluding pieces not included in this part
        % by inserting NaNs between intervals we can plot this as one line, which saves on memory
        dindax = abs([0; diff(ppot)])>0.1;
        pos_plot = [ppox ppoy];
        pos_plot(dindax,:) = NaN;
        plot(pos_plot(:,1),pos_plot(:,2),'Color',[.5 .5 .5 .5]); hold on;

        % plot spikes after position data so they are all on top     
        plot(pspx,pspy,'Color',[1 0 0 0.5],'Marker','.','LineStyle','none') 
        
        % additional settings
        daspect([1 1 1])
        axis xy off tight
        text(0,1.1,sprintf('%d spikes (%.2f Hz), %d seconds (%.1f mins)',numel(pspx),numel(pspx)/part_duration,round(part_duration),part_duration/60),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
        
%% >>>>>>>>>> Firing rate map       
    axrt = axes('Units','pixels','Position',[xvec(cc),yvec(2),pwidth,pheight]);
        im = imagesc(ratemap,'alphadata',~isnan(ratemap));
        daspect([1 1 1])
        caxis([0 max([0.1 max(ratemap(:),[],'omitnan')])])       
        %colormap(axrt,flipud(viridis));
        colormap(axrt,turbo)
        axis xy off
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);   
 
        text(0,1.05,sprintf('SI: %.2f, z: %.2f',sinfo(1),sinfo(2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')

        axp = get(gca,'Position');
        ccb = colorbar;
        set(gca,'Position',axp);
        set(ccb,'Position',get(ccb,'Position')+[0 0 -0.004 0])
        title(ccb,'Hz','FontSize',fsiz)        
        set(ccb,'yticklabel',num2str(get(ccb,'ytick')','%.1f')) % change Ytick labels to have the same number of decimal places
                
%% >>>>>>>>>> Grid autocorrelogram              
    axgs = axes('Units','pixels','Position',[xvec(cc),yvec(3),pwidth,pheight]);
        imc = imagesc(amap,'alphadata',~isnan(amap));
        daspect([1 1 1])
        caxis([-0.2 1])
        colormap(axgs,turbo);        
        axis xy off
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
        
        text(0,1.05,sprintf('G: %.2f, z: %.2f',sinfo(3),sinfo(4)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')

        axp = get(gca,'Position');
        ccb = colorbar;
        set(gca,'Position',axp);
        set(ccb,'Position',get(ccb,'Position')+[0 0 -0.004 0])
        title(ccb,'Hz','FontSize',fsiz)        
        set(ccb,'yticklabel',num2str(get(ccb,'ytick')','%.1f')) % change Ytick labels to have the same number of decimal places
                       
    end

%% >>>>>>>>>> Save the overall figure
        % Save the figure  
        % [~,~,~] = mkdir([pdata.directory '\' pdata.outname '\part_figures']); % create a folder to hold outputs 
        % fname = [pdata.directory '\' pdata.outname '\part_figures\' uci '_projections.png'];
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




































