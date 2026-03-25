function add_manual_cell_type(opts)
% name short description
% longer description
%
% USAGE
%
% out = name(in) process with default settings
%
% out = name(in,optional) process using optional argument 1
%
% out = rate_mapper(__,name,value) process with Name-Value pairs 
%
% INPUT
%
% 'in' - Scalar, positive integer that specifies X
%       units are in Y.
%       Default value is Z.
%
% OUTPUT
%
% 'out' - Scalar, positive integer that specifies X
%       units are in Y.
%       Default value is Z.
%
% NOTES
% 1. 
%
% 2. 
%
% EXAMPLE
% 
% SEE ALSO NAME, NAME

% HISTORY
%
% version 1.0.0, Release 00/00/26 Initial release
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
%%%%%%%%%%%%%%%% ARGUMENT CHECK
    arguments
        opts.combine_sdata          logical = false
        opts.master_sdata           double = []
        opts.master_bdata           double = []  
        opts.data_dir               string = pwd
        opts.klustest               logical = false
        opts.template               logical = false
    end

    % settings
    figvis = 'on';
    global cell_type;         
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    fig_overall = figure('visible',figvis,'Units','pixels','Position',[20, 150, 1800, 700]); % open/create figure
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    fsiz = 8; % the fontsize for different texts
    ann = annotation('textbox',[0, 1, 1, 0],'string','','FontSize',8,'LineStyle','none','interpreter','none');  

%%%%%%%%%%%%%%%% Prepare the data
    sdata_dir = fullfile(opts.data_dir,'klustest','sdata.mat');
    load(sdata_dir,'sdata');

    if ~any(ismember(sdata.Properties.VariableNames,'cell_type')) % if the column(s) do not exist yet
        sdata.cell_type = NaN(size(sdata,1),2); % preallocate
    end

    pdata = sdata.Properties.CustomProperties.pdata;
    part_names = fieldnames(pdata.part_config);
    bdata = sdata.Properties.CustomProperties.bdata;

    % Run through every cluster
    ucis = unique(sdata.uci); % list of unique cells in sdata  
    for uu = 1:length(ucis)
        uci = ucis{uu};
        disp(sprintf('\tCell %d of %d (%.f%%): %s',uu,length(ucis),uu/length(ucis)*100,uci))
        
        % clear figure for next cluster
        clf(fig_overall);     
        delete(ann);
        
        % add an annotation to the figure with some important info
        ann_str = sprintf('Cell: %s, Analysed: %s',uci,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        ann = annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',12,'LineStyle','none','interpreter','none');  
        
        idx = find(ismember(sdata.uci,uci));
        
        % get info for automated cell typing
        [~,midx] = max(sdata.frate(idx));
        idx = idx(midx);
        
        frate = sdata.frate(idx); % firing rate
        wow = sdata.wave_width(idx); % width of waveform, this is for the channel with the highest mean amplitue
        spat = sdata.spatial_info(idx,1); % spatial information
        gscore = sdata.grid_info(idx,1); % grid score
        rayv = sdata.hd_info(idx,1); % rayleigh vector
        rpv = sdata.autocorr_refrac_info(idx,1);
        rpvprop = sdata.autocorr_refrac_info(idx,2);

        % do automated cell typing
        rate_min = 0.1;
        rate_max = 10;
        width_cutoff = 0.250;
        spatial_cutoff = 1; % b/s
        grid_cutoff = 0.8; % g-score
        hd_cutoff = 0.2; % rayleigh
        rpv_cutoff = 5; % percentage cutoff for refractory period violations
        % 1 = noise
        % 2 = pyramidal
        % 3 = place cell
        % 4 = grid cell
        % 5 = hd cell
        % 6 = interneuron
        cell_type = 1; % noise by default
        cols = {'g','w','w','w','w','w','w'}; % for plotting later        
        if frate > rate_min && rpvprop < rpv_cutoff % if the cluster's firing is not too low and there are not too many RPVs
            if frate < rate_max % if the cell's firing rate does not exceed the maximum
                if wow > width_cutoff % if the cell has a pyramidal waveform
                    cell_type = 2; % pyramidal
                    cols = {'w','w','w','w','w','w','w'}; % for plotting later        
                    cols{2} = 'g';
                    if spat>spatial_cutoff % if the pyramidal cell is spatial
                        cell_type = 3; % place cell
                        cols = {'w','w','w','w','w','w','w'}; % for plotting later                              
                        cols{3} = 'g';
                    end                     
                else
                    cell_type = 6; % low firing interneuron
                    cols = {'w','w','w','w','w','w','w'}; % for plotting later                          
                    cols{6} = 'g';
                end   
            else
                cell_type = 6; % high firing interneuron
                cols = {'w','w','w','w','w','w','w'}; % for plotting later                          
                cols{6} = 'g';                
            end
            if gscore>grid_cutoff % if the grid score is significant
                cell_type = 4; % high firing interneuron     
                cols = {'w','w','w','w','w','w','w'}; % for plotting later                          
                cols{4} = 'g';                
            end 
            if rayv>hd_cutoff % if the rayleigh vector is significant
                cell_type = 5; % high firing interneuron     
                cols = {'w','w','w','w','w','w','w'}; % for plotting later                          
                cols{5} = 'g';                
            end
        end
        cval1 = cell_type; % automated value
        
%%%%%%%%%%%%%%%% Manual refining
        % pushbuttons for cell typing
        bxpos = 1690;
        ymx = 600;
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{1},'String','Noise','Position',[bxpos ymx 100 50],'Callback',{@change_ctype,1});  
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{2},'String','Pyramidal cell','Position',[bxpos ymx-80 100 50],'Callback',{@change_ctype,2}); 
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{3},'String','Place cell','Position',[bxpos ymx-160 100 50],'Callback',{@change_ctype,3});   
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{7},'String','Spatial cell','Position',[bxpos ymx-240 100 50],'Callback',{@change_ctype,7});           
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{6},'String','Interneuron','Position',[bxpos ymx-320 100 50],'Callback',{@change_ctype,6});              
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor','m','String','Agree','Position',[bxpos ymx-400 100 50],'Callback',{@change_ctype,cell_type});              
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{4},'String','Grid cell','Position',[bxpos ymx-480 100 50],'Callback',{@change_ctype,4});     
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{5},'String','HD cell','Position',[bxpos ymx-560 100 50],'Callback',{@change_ctype,5});  
        
        xmin = 20;
        ymax = 450;
        ybuff = 200;
        yvec = [ymax ymax-ybuff ymax-2*ybuff ymax-3*ybuff ymax-4*ybuff ymax-5*ybuff];
        pwidth = 220;
        pheight = 190;        
        
        nparts = length(part_names);
        for pp = 1:nparts
            idx = find(ismember(sdata.uci,uci) & sdata.partn==pp);
            if isempty(idx)
                continue
            end
            if isempty(sdata.spt_pot_index{idx})
                continue
            end

            % part_duration = sdata.session_duration_s(idx);
            part_now = part_names{pp};
            pos = pdata.pos;
            pidx = pos.pot > pdata.part_config.(part_now).interval_times(:,1) & pos.pot < pdata.part_config.(part_now).interval_times(:,2); % index for position data in this part
            ppox = pos.pox(pidx);
            ppoy = pos.poy(pidx);
            ppot = pos.pot(pidx);
            pspx = pos.pox(sdata.spt_pot_index{idx});
            pspy = pos.poy(sdata.spt_pot_index{idx});
            ratemap = sdata.ratemap{idx};
            amap = sdata.gridmap{idx};           

%%%%%%%%%%%%%%%% Positions and spikes
            xnow = xmin;
            axps = axes('Units','pixels','Position',[xnow,yvec(pp),pwidth,pheight]);
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
                part_duration = numel(ppot).*pdata.pos_srate;
                text(0,1.1,sprintf('%d spikes (%.2f Hz), %.1f mins',numel(pspx),numel(pspx)/part_duration,part_duration/60),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
        
%%%%%%%%%%%%%%%% Firing rate map       
            xnow = xnow + pwidth + 20;
            axrt = axes('Units','pixels','Position',[xnow,yvec(pp),pwidth,pheight]);
                im = imagesc(ratemap,'alphadata',~isnan(ratemap));
                daspect([1 1 1])
                caxis([0 max([0.1 max(ratemap(:),[],'omitnan')])])       
                colormap(axrt,turbo);
                axis xy off

                sinfo = sdata.spatial_info(idx,:);            
                text(0,1.05,sprintf('SI: %.2f',sinfo(1)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')

                axp = get(gca,'Position');
                ccb = colorbar;
                set(gca,'Position',axp);
                set(ccb,'Position',get(ccb,'Position')+[0 0 -0.004 0])
                title(ccb,'Hz','FontSize',fsiz)        
                set(ccb,'yticklabel',num2str(get(ccb,'ytick')','%.1f')) % change Ytick labels to have the same number of decimal places
                
%%%%%%%%%%%%%%%% Grid autocorrelogram     
            xnow = xnow + pwidth + 20;
            axgs = axes('Units','pixels','Position',[xnow,yvec(pp),pwidth,pheight]);
                imc = imagesc(amap,'alphadata',~isnan(amap));
                daspect([1 1 1])
                caxis([-0.2 1])
                colormap(axgs,turbo);        
                axis xy off

                ginfo = sdata.grid_info(idx,:);            
                text(0,1.1,sprintf('G: %.2f',ginfo(1)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')

%%%%%%%%%%%%%%%% Head direction  
            xnow = xnow + pwidth + 20;
            axhd = axes('Units','pixels','Position',[xnow,yvec(pp),150,180]);
                hd_ratemap = sdata.hd_ratemap{idx};
                hd_dwellmap = pdata.(part_now).hd_dwellmap;

                theta = linspace(0,2*pi,pdata.mapset.hd_bins)'; 
                p1 = polarplot(theta,hd_ratemap,'k'); hold on;
                dwellmap_to_plot = (hd_dwellmap/max(hd_dwellmap(:)))*max(hd_ratemap)*0.5;
                p1 = polarplot(theta,dwellmap_to_plot,'b'); hold on;            
                axhd2.ThetaZeroLocation = 'right';
    
                axhd2.Color = 'white';
                axhd2.GridColor = 'black'; % Change grid line color to black
                axhd2.ThetaColor = 'black'; % Change theta axis color to black
                axhd2.RColor = 'black'; % Change r axis color to black
    
                text(0,1.2,sprintf('r: %.2f, PFD%c: %.2f',sdata.hd_info(idx,1),176,sdata.hd_info(idx,2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')

%%%%%%%%%%%%%%%% Waveforms      
            xnow = xnow + pwidth + 20;
            axwav = axes('Units','pixels','Position',[xnow,yvec(pp),180,180]);
                mxs = max(sdata.wave_mean{idx},[],2,'omitnan');
                [~,widx] = sort(mxs,'descend'); % sort from largest > smallest waveform
                wavtime = pdata.wavtime; 

                [hl,hp] = boundedline(wavtime,sdata.wave_mean{idx}(widx(1),:),sdata.wave_std{idx}(widx(1),:),'-k');
                set(hl,'Color','r','LineStyle','-','LineWidth',1); 
                set(hp,'FaceColor','b','FaceAlpha',.5);     
                ax = gca;
                ax.XLim = [-0.2 0.8] .* 1e3;
                ax.YDir = 'normal';
                box on
                grid on
                if pp<3
                    ax.XTick = [];  
                else
                    xlabel('Time (ms)')
                end                
                ylabel(sprintf('Amplitude (%cV)',char(956)))
                text(1,1.1,sprintf('Ch%d, Peak: %.1f%cV, Width: %.2f%cs',widx(1),mxs(widx(1)),char(956),sdata.wave_width(idx).*1000,char(956)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','right')

%%%%%%%%%%%%%%%% Spike autocorrelation - 25ms refractory period   
            xnow = xnow + pwidth + 20;
            axref = axes('Units','pixels','Position',[xnow,yvec(pp),180,180]); 
                f = sdata.autocorr_refrac{idx};
                xi = pdata.autocorr_refrac_xvalues;
    
                edg = pdata.autocorr_refrac_evalues;
                x = xi(:)';
                ex = edg(:)';
                y = f(:)';
    
                % marker XY matrices
                X = [ex(1:end-1); ex(2:end); ex(2:end); ex(1:end-1)];
                Y = [zeros(size(y)); zeros(size(y)); y; y];   
                C = zeros(numel(y),3);
                C(abs(x(:))<2,:) = ones(sum(abs(x(:))<2),3).*[1 0 0];
            
                % plot the data 
                patch(axref,X,Y,'k','facevertexcdata',C,'facecolor','flat','edgecolor','none'); % plot from scratch
    
                % axis settings
                ax = gca;            
                ax.XLim = [-20 20];
                ax.XTick = -20:10:20;            
                xlabel('Time lag (ms)') % label x axis
                ylabel('Probability') % label y axis
                set(gca,'LineWidth',1,'layer','top','FontSize',fsiz); 

                text(0,1.05,sprintf('RPV: %.2f, RPV%%: %.2f',rpv,rpvprop),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
                set(ax,'yticklabel',num2str(get(gca,'ytick')','%.3f'))        
                grid on;

%%%%%%%%%%%%%%%% Spike autocorrelation - 500ms theta modulation     
            xnow = xnow + pwidth + 20;
            axref = axes('Units','pixels','Position',[xnow,yvec(pp),200,180]);
                f = sdata.autocorr_theta{idx}(:,1);
                xi = pdata.autocorr_theta_xvalues;
    
                edg = pdata.autocorr_theta_evalues;
                x = xi(:)';
                ex = edg(:)';
                y = f(:)';
    
                % marker XY matrices
                X = [ex(1:end-1); ex(2:end); ex(2:end); ex(1:end-1)];
                Y = [zeros(size(y)); zeros(size(y)); y; y];   
                C = zeros(numel(y),3);
                C(abs(x(:))<2,:) = ones(sum(abs(x(:))<2),3).*[1 0 0];
            
                % plot the data 
                patch(axref,X,Y,'k','facevertexcdata',C,'facecolor','flat','edgecolor','none'); % plot from scratch
    
                % axis settings
                ax = gca;
                ax.XLim = [-500 500];
                ax.XTick = -500:125:500;
                ylabel('Probability') % label y axis
                set(gca,'LineWidth',1,'layer','top','FontSize',fsiz); 
                set(ax,'yticklabel',num2str(get(gca,'ytick')','%.3f'))        
                grid on;          
                xlabel('Time lag (ms)') % label x axis        

        end

%%%%%%%%%%%%%%%% get users response and add it to sdata
        uiwait(fig_overall); % at this point start waiting for a cell type button to be pressed

        cval2 = cell_type; % curated value
        sdata.cell_type(ismember(sdata.uci,uci),:) = repmat([cval1 cval2],sum(ismember(sdata.uci,uci)),1);  % automated cell type, manual cell type 

    end
    close(fig_overall)

%%%%%%%%%%%%%%%% %% Save the data and finish up
    save(sdata_dir,'sdata'); % save session data
    analysis_log({'add_manual_cell_type'},1,'version',{'v1.0.0'});
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pushbutton control function
function change_ctype(~,~,value) 
    global cell_type; 
    cell_type = value; 
    uiresume(gcbf);
end









































