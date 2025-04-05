function kluspartfig(sdata,uci,fig_vis,save_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an sdata structure, a cell identifier and generates a figure with plots for each part identified
%   kluspartfig(sdata,uci,part_now,fig_vis,save_fig)
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
if ~exist('fig_vis','var') || isempty(fig_vis)
    fig_vis = 'on';
end % if ~exist('fig_vis','var') || ismepty(fig_vis)

if ~exist('save_fig','var') || isempty(save_fig)
    save_fig = 0;
end % if ~exist('save_fig','var') || ismepty(save_fig)

part_names = cell(1,length(sdata.part_config));
for pp = 1:length(sdata.part_config)
    pname = sdata.part_config{pp}{1};
    part_names{pp} = pname;
end % for pp = 1:length(sdata.part_config)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create figure
fig_overall = figure('visible',fig_vis,'Position',[100, 100, 1200, 800]);
set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white
colormap(jet(256)); % to make sure the colormap is not the horrible default one
fig_hor = numel(part_names)+1; % how many plots wide should it be
fig_ver = 3; % how many plots tall should it be
fspac = 0.01; % the spacing around the plots, on all sides
fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
fmarg = 0.03; % the margins around the plots, at the edge of the figure
fsiz = 8; % the fontsize for different texts
flw = 1; % the line width for different plots

% add an annotation to the figure with some important info
ann_str = sprintf('Cell: %s, Analysed: %s',uci,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');  

maps = cell(numel(part_names));
for pp = 1:numel(part_names)
    part_now = part_names{pp};
    part_duration = sdata.(part_now).duration; 
    pox = sdata.(part_now).pox;
    poy = sdata.(part_now).poy;
    pot = sdata.(part_now).pot;
    spx = sdata.(uci).(part_now).spx;
    spy = sdata.(uci).(part_now).spy;
    ratemap = sdata.(uci).(part_now).ratemap;
    skaggs = sdata.(uci).(part_now).spatial_measures.spatial_information;
    spars = sdata.(uci).(part_now).spatial_measures.sparsity;
    MMF = sdata.(uci).(part_now).spatial_measures.mean_method_focus;

    %% spike plot with black lines for path and red dots for spikes
    subaxis(fig_ver,fig_hor,pp,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
    plot(pox,poy,'k')
    hold on
    plot(spx,spy,'ro','MarkerFaceColor','r','MarkerSize',2)
    daspect([1 1 1])
    axis xy off
    axis([min(pox)-15 max(pox)+15 min(poy)-15 max(poy)+15]);
    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
    title(sprintf('%s %d spikes (%.2f Hz)',part_now,numel(spx),numel(spx)/part_duration),'FontSize',6);

    %% firing rate map
    subaxis(fig_ver,fig_hor,pp+numel(part_names)+1,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
    im = imagesc(ratemap);
    set(im,'alphadata',~isnan(ratemap));
    daspect([1 1 1])
    axis xy off
    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);   
    title(sprintf('Skaggs: %.2f, Sparsity: %.2f, MMF: %.2f',skaggs,(spars*100),MMF),'FontSize',6);
    tt = text(-1,5,sprintf('Mx: %.2f (Hz), Mu: %.2f (Hz), Mn: %.2f (Hz)',nanmax(ratemap(:)),nanmean(ratemap(:)),nanmin(ratemap(:))),'FontSize',6);
    set(tt,'rotation',90);
    maps{pp} = ratemap;
    
end % for pp = 1:numel(part_names)

%% map correlation plot
epairs = nchoosek(1:numel(part_names),2); % every possible combination of ratemap
epairs = [epairs; fliplr(epairs)];
cmat = NaN(numel(part_names),numel(part_names));
for cc = 1:length(epairs(:,1))
    m1 = maps{epairs(cc,1)};
    m2 = maps{epairs(cc,2)};
    if size(m1) ~= size(m2)
        m2 = imresize(m2,size(m1));
    end % if size(m1) ~= size(m2)
    [r,p] = corr(m1(:),m2(:),'Type','Spearman','rows','pairwise');
    cmat(epairs(cc,1),epairs(cc,2)) = r;
end % for cc = 1:length(epairs(:,1))

subaxis(fig_ver,fig_hor,numel(part_names)+1,'Spacing',fspac,'Padding',fpadd,'Margin',0.05,'Holdaxis');
im = imagesc(cmat);
set(im,'alphadata',~isnan(cmat));
axis square on
ax = gca;
ax.XTick = 1:numel(part_names);
ax.YTick = 1:numel(part_names);
ax.XTickLabel = part_names;
ax.YTickLabel = part_names;
ax.YTickLabelRotation = 90;

xlabel('Part 1');
ylabel('Part 2');
title('Resized map correlations (Spearman pairwise)','FontSize',6);
for cc = 1:length(epairs(:,1))
    text(epairs(cc,2),epairs(cc,1),sprintf('%.1f',cmat(epairs(cc,1),epairs(cc,2))),'Color','w','HorizontalAlignment','center','VerticalAlignment','middle')
end % for cc = 1:length(epairs(:,1))

%% Spikes vs time plot
axf = subaxis(fig_ver,fig_hor,(numel(part_names)*2+3):(numel(part_names)*3+3),'Spacing',fspac,'Padding',fpadd,'Margin',0.05,'Holdaxis');
colormap(axf,jet);
pot = sdata.(uci).spike_times;
bstvals = (min(pot):sdata.config.time_bins:max(pot)); % vector of time points at which we should calculate spike probability
[bspikes,~] = histc(pot,bstvals);
[bsprobs,xibs] = ksdensity(pot,bstvals);

yyaxis left
bar(bstvals,bspikes,1,'FaceColor','k');
ylabel('Spikes')
v = axis;
hold on
ax = gca;
ax.YColor = 'k';

cols = copper(numel(part_names));
for pp = 1:numel(part_names)
    part_now = part_names{pp};
    potn = sdata.(part_now).pot;
    patch([min(potn) min(potn) max(potn) max(potn)],[v(3) v(4) v(4) v(3)],cols(pp,:),'FaceAlpha',0.2,'EdgeColor','none')
    text(mean([max(potn) min(potn)]),v(4)-0.05*v(4),part_now)
end % for pp = 1:numel(part_names)
set(ax,'children',flipud(get(gca,'children')))

yyaxis right
plot(xibs,bsprobs,'r');
ylabel('Probability')

ax = gca;
ax.XLim = [min(pot) max(pot)];
ax.YColor = 'r';
xlabel('Time (s)')











