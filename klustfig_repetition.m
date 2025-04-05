function klustfig_repetition(pdata,cluma,rnow,dnow,pp,fast_figs,figvis,fname,part_now)
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
        ann_str = sprintf('Rat: %s, Date: %s, Part: %s, Analysed: %s',rnow,dnow,part_now,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',8,'LineStyle','none','interpreter','none');  
            
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> ranked autocorrelation centers
        ax = axes('Units','pixels','Position',[50,100,200,650]);
            pidx = cluma.partn==pp;
            if any(ismember(cluma.Properties.VariableNames,'cell_type')) % if this session has been cell typed
                ctype = cluma.cell_type(pidx,1); % load them
                pidx(ctype~=3) = false; % remove anything that is not a place cell
            end  
            dat = cat(1,cluma.amap_cent{pidx});
            amaps = cat(3,cluma.planar_amap{pidx}); 
            scores = cluma.repetition_score(pidx,1);            
            [~,idx] = sort(scores,'descend','MissingPlacement','last');
            dat = dat(idx,:);
            scores = scores(idx);
            amaps = amaps(:,:,idx);
            row_idx = find(pidx);
            row_idx = row_idx(idx);

            imagesc(flipud(dat))
            colormap(gca,'turbo')
            axis xy
            xlabel('Lag (bins)')
            ylabel('Ranked cell')
            caxis([-0.2 1])
    
%% >>>>>>>>>> Top N cells         
        xmin = 300;
        ymax = 650;
        xsiz = 180;
        ysiz = 120;
        xbuff = 15;
        ybuff = 0;
        xplots = 6;
        yplots = 6;
        xvec = xmin+(0:xplots-1)+((0:xplots-1)*xsiz)+((0:xplots-1)*xbuff);
        yvec = ymax-(0:yplots-1)-((0:yplots-1)*ysiz)-((0:yplots-1)*ybuff);
        [xx,yy] = meshgrid(xvec,yvec);
        xx = xx;
        yy = yy;
        
        for ff = 1:min([sum(pidx) xplots*yplots])
            ax = axes('Units','pixels','Position',[xx(ff),yy(ff),xsiz,ysiz]);
                map = cluma.ratemap_planar{row_idx(ff)};
                imagesc(map,'alphadata',~isnan(map));
                axis xy
                daspect([1 1 1])
                %colormap(ax,flipud(viridis))
                colormap(gca,'turbo')   
                caxis( [0 max([max(map(:),[],'omitnan') 1])] )                
                axis off
                
                text(0,1.15,sprintf('Score = %.2f\n%s',scores(ff),cluma.uci{row_idx(ff)}),'Units','normalized','HorizontalAl','left','Interpreter','none','FontSize',8)
        end
        
%% >>>>>>>>>> Save the overall figure
    % Save the figure  
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
    end

 




































