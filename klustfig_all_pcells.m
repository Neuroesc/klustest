function klustfig_all_pcells(pdata,cluma,rnow,dnow,pp,fast_figs,figvis,fieldz,fname)
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
        fig_clust = figure('visible',figvis,'Units','pixels','Position',[10, 50, 1810, 900]); % open/create figure
        set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
        set(gcf,'color','w'); % makes the background colour white
        colormap(jet(256)); % to make sure the colormap is not the horrible default one
        fsiz = 8; % the fontsize for different texts
        flnw = 0.5; % the line width for different plots

        % add an annotation to the figure with some important info
        part_names = {'arena1','hills','arena2'};        
        part_now = part_names{pp};        
        ann_str = sprintf('Rat: %s, Date: %s, Part: %s, Analysed: %s',rnow,dnow,part_now,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',8,'LineStyle','none','interpreter','none');  
            
        if ~exist('fieldz','var')
            fieldz = 0; % 1 = plot place fields mode
        end
        
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> ranked autocorrelation centers
        idx = cluma.session==pp;
        if any(ismember(cluma.Properties.VariableNames,'cell_type')) % if this session has been cell typed
            ctype = cluma.cell_type(idx,1); % load them
            idx(ctype~=3) = false; % remove anything that is not a place cell
        end  
        scores = cluma.planar_spatial_info_shuffles(idx,1);        
        [~,idxn] = sort(scores,'descend','MissingPlacement','last');
        row_idx = find(idx);
        row_idx = row_idx(idxn);
        
%% >>>>>>>>>> Top N cells         
        xmin = 50;
        ymax = 770;
        xsiz = 180;
        ysiz = 90;
        xbuff = 10;
        ybuff = 4;
        xplots = 9;
        yplots = 9;
        xvec = xmin+(0:xplots-1)+((0:xplots-1)*xsiz)+((0:xplots-1)*xbuff);
        yvec = ymax-(0:yplots-1)-((0:yplots-1)*ysiz)-((0:yplots-1)*ybuff);
        [xx,yy] = meshgrid(xvec,yvec);
        
        for ff = 1:min([numel(row_idx) xplots*yplots])
            ax = axes('Units','pixels','Position',[xx(ff),yy(ff),xsiz,ysiz]);            
                if fieldz
                    map = cluma.ratemap_planar{row_idx(ff)};
                    mapset = pdata.mapset;
                    
                    zmap = ( map - mean(map(:),'omitnan') ) / std(map(:),'omitnan'); % zscore ratemap
                    thresh_ratemap = imbinarize(zmap,mapset.zcut); % 2 s.d. threshold
                    datout = regionprops('table',thresh_ratemap,map,'Area','PixelIdxList','MaxIntensity');
                    bmap = ones(size(map)).*0.25;                                        
                    if ~isempty(datout)
                        datout.Area(:) = datout.Area(:) .* ((mapset.binsize/10)^2); % convert field area to cm2
                        nindx = datout.Area(:) < mapset.arcut | datout.MaxIntensity(:) < mapset.frcut;                        
                        datout(nindx,:) = [];  
                        for kk = 1:size(datout,1)
                            pixidx = datout.PixelIdxList{kk};
                            bmap(pixidx) = 1;
                        end                        
                    end                    
                    imagesc(map,'alphadata',bmap);
                    axis xy
                    daspect([1 1 1])
                    colormap(ax,turbo)   
                    caxis( [0 max([max(map(:),[],'omitnan') 1])] )
                    axis off
                    text(0,1,sprintf('t%dc%d %.1fb/s %dfields',cluma.ele(row_idx(ff)),cluma.clu(row_idx(ff)),cluma.planar_spatial_info_shuffles(row_idx(ff),1),size(datout,1)),'Color','k','Units','normalized','HorizontalAl','left','Interpreter','none','FontSize',8,'VerticalAl','top')

                else
                    map = cluma.ratemap_planar{row_idx(ff)};
                    imagesc(map,'alphadata',~isnan(map));
                    axis xy
                    daspect([1 1 1])
                    %colormap(ax,flipud(viridis))
                    colormap(ax,turbo)   
                    caxis( [0 max([max(map(:),[],'omitnan') 1])] )
                    axis off
                    text(0,1,sprintf('t%dc%d %.1fHz',cluma.ele(row_idx(ff)),cluma.clu(row_idx(ff)),max(map(:),[],'omitnan')),'Color','w','Units','normalized','HorizontalAl','left','Interpreter','none','FontSize',8,'VerticalAl','top')
                end
        end

%% >>>>>>>>>> Save the overall figure
    % Save the figure  
    % [~,~,~] = mkdir([pdata.directory '\' pdata.outname]); % create a folder to hold outputs 
    % if fieldz
    %     fname = [pdata.directory '\' pdata.outname '\' part_now '_all_place_cell_fields.png'];        
    % else
    %     fname = [pdata.directory '\' pdata.outname '\' part_now '_all_place_cells.png'];
    % end
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


 




































