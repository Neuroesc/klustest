function figPARTS(sdata,uci,fig_vis,save_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an sdata structure, a cell identifier and generates a figure with plots for each part identified
%   figPARTS(sdata,uci,part_now,fig_vis,save_fig)
%
%%%%%%%% Inputs
%   sdata = sdata structure
%   uci = unique cell identifier
%   fig_vis = (optional), 'on' to show figures, 'off' to hide figures, default is 'off'
%   save_fig = (optional), 1 to save figures in .fig files, 0 otherwise, default is 0
%
%%%%%%%% Comments
%   02/04/17 created to contain all of this figure code
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
close all

if ~exist('fig_vis','var') || isempty(fig_vis)
    fig_vis = 'on';
end % if ~exist('fig_vis','var') || ismepty(fig_vis)

if ~exist('save_fig','var') || isempty(save_fig)
    save_fig = 0;
end % if ~exist('save_fig','var') || ismepty(save_fig)

part_names = fieldnames(sdata.part_config);
part_names = part_names(1:end-1);

% determine cluster and tetrode
sst = strsplit(uci,'_');
tet = str2double(sst{3}(regexp(sst{3},'\d')));
tetstr = ['t' num2str(tet)];
clu = str2double(sst{4}(regexp(sst{4},'\d')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create figure
fig_overall = figure('visible',fig_vis,'Position',[100, 100, 1600, 800]);
set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white
colormap(jet(256)); % to make sure the colormap is not the horrible default one
fig_hor = numel(part_names); % how many plots wide should it be
fig_ver = 3; % how many plots tall should it be
fspac = 0.01; % the spacing around the plots, on all sides
fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
fmarg = 0.03; % the margins around the plots, at the edge of the figure
fsiz = 8; % the fontsize for different texts
flnw = 1; % the line width for different plots

% add an annotation to the figure with some important info
ann_str = sprintf('Cell: %s, Analysed: %s',uci,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');  

maps = cell(numel(part_names));
for pp = 1:numel(part_names)
    part_now = part_names{pp};
    part_duration = sdata.(part_now).duration; 
    
    % spatial data
    pox = sdata.(part_now).pox; % position x for this part
    poy = sdata.(part_now).poy; % position y for this part
    pot = sdata.(part_now).pot; % position time for this part
    spx = sdata.(uci).(part_now).spx; % spike x for this part
    spy = sdata.(uci).(part_now).spy; % spike y for this part
    spt = sdata.(uci).(part_now).spt; % spike t for this part  
    
    if ~numel(spx) % if there are no spikes
        maps{pp} = NaN;
        continue
    end
    
    ratemap = sdata.(uci).(part_now).ratemap;
    skaggs = sdata.(uci).(part_now).spatial_measures.spatial_information;
    spars = sdata.(uci).(part_now).spatial_measures.sparsity;
    MMF = sdata.(uci).(part_now).spatial_measures.mean_method_focus;

    %% spike plot with black lines for path and red dots for spikes
    subaxis(fig_ver,fig_hor,pp,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
    tvals = sdata.part_config.(part_now).times;
    hold on
    for tt = 1:length(tvals(:,1))
        tindx = pot > tvals(tt,1) & pot < tvals(tt,2);
        plot(pox(tindx),poy(tindx),'k')
    end    

    % plot spikes after so they are always on top
    for tt = 1:length(tvals(:,1))
        sindx = spt > tvals(tt,1) & spt < tvals(tt,2);
        plot(spx(sindx),spy(sindx),'ro','MarkerFaceColor','r','MarkerSize',2)
    end   
    daspect([1 1 1])
    axis xy off tight square
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
    title(sprintf('%s %d spikes (%.2f Hz)',part_now,numel(spx),numel(spx)/part_duration),'FontSize',6);

    %% firing rate map
    subaxis(fig_ver,fig_hor,pp+numel(part_names),'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
    im = imagesc(ratemap);
    set(im,'alphadata',~isnan(ratemap));
    daspect([1 1 1])
    axis xy off tight square
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);   
    title(sprintf('Skaggs: %.2f, Sparsity: %.2f, MMF: %.2f',skaggs,(spars*100),MMF),'FontSize',6);
    tt = text(-1,5,sprintf('Mx: %.2f (Hz), Mu: %.2f (Hz), Mn: %.2f (Hz)',nanmax(ratemap(:)),nanmean(ratemap(:)),nanmin(ratemap(:))),'FontSize',6);
    set(tt,'rotation',90);
    maps{pp} = ratemap;
    
end 

%% Spikes vs time plot
axst = axes('Units','pixels','Position',[50,120,1500,80]);
    spt = double(sdata.(uci).spike_times);
    pot = double(sdata.pot); 
    pov = double(sdata.pov); 

    plot(pot,pov,'k')
    ylabel('Speed (cm/s)')
    ax = gca;
    ax.XTick = [];
    ax.XLim = [min(pot) max(pot)];

    cols = copper(numel(part_names));
    for pp = 1:numel(part_names)
        part_now = part_names{pp};
        potn = double(sdata.(part_now).pot);
        patch([min(potn) min(potn) max(potn) max(potn)],[ax.YLim ax.YLim(2:-1:1)],cols(pp,:),'FaceAlpha',0.2,'EdgeColor','none')
        text(mean([max(potn) min(potn)]),ax.YLim(2)-0.05*ax.YLim(2),part_now)
    end
    set(ax,'children',flipud(get(gca,'children')))
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % set the axes to be on top, set their thickness to flnw, set all fonts to size fsiz and make sure the axis and text are black

axst = axes('Units','pixels','Position',[50,35,1500,80]);
    bstvals = 0:sdata.config.time_bins:sdata.session_duration; % vector of time points at which we should calculate spike probability
    cents = histcents(bstvals);
    [bspikes,~] = histcounts(spt,bstvals);
    [bsprobs,xibs] = ksdensity(spt,bstvals);

    yyaxis left
    bar(cents,bspikes,1,'FaceColor','k');

    ylabel('Spikes')
    hold on
    ax = gca;
    ax.YColor = 'k';
    ax.YLim(2) = nanmax(bspikes)+1;

    cols = copper(numel(part_names));
    for pp = 1:numel(part_names)
        part_now = part_names{pp};
        potn = double(sdata.(part_now).pot);
        patch([min(potn) min(potn) max(potn) max(potn)],[ax.YLim ax.YLim(2:-1:1)],cols(pp,:),'FaceAlpha',0.2,'EdgeColor','none')
    end
    set(ax,'children',flipud(get(gca,'children')))

    yyaxis right
    plot(xibs,bsprobs,'r');
    ylabel('Probability')

    ax = gca;
    ax.XLim = [0 sdata.session_duration];
    ax.YColor = 'r';
    xlabel('Time (s)')
    set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % set the axes to be on top, set their thickness to flnw, set all fonts to size fsiz and make sure the axis and text are black

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the figure    
figfile = [pwd '\klustest\' sdata.settings.combined_name '\figures\'];
[~,~,~] = mkdir(figfile);
print(fig_overall,'-dpng','-r150',[figfile uci '_parts.png'])
if save_fig
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')');
    savefig(fig_overall,[figfile uci '_parts.fig'],'compact');
end
close(fig_overall);                                          
















