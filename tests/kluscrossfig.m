function kluscrossfig(tet,clusters,sdata,fig_vis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function generates a figure showin the cross- spike correlogram for clusters on a tetrode
%   kluscrossfig(tet,clusters,mtint,sdata,fig_vis)
%
%%%%%%%% Inputs
%   tet = the tetrode to look at
%   clusters = the cluster vector on this tetrode
%   mtint = an mtint structure
%   sdata = an sdata structure
%   fig_vis = figure visibility, set to 'on' to see figure, 'off' to hide figure
%
%%%%%%%% Comments
%   31/03/17 created to contain all of this figure code
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('fig_vis','var') || isempty(fig_vis)
    fig_vis = 'off';
end % if ~exist('fig_vis','var') || ismepty(fig_vis)

clusters(clusters == 0) = [];
clus_count = numel(clusters);
spik_count = sdata.spike_count{tet};
duration = sdata.session_duration;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot figure
fig_corr = figure('visible',fig_vis,'Position',[100, 100, 1024, 800]);
set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white
colormap(jet(256)); % to make sure the colormap is not the horrible default one
fig_hor = clus_count+1; % how many plots wide should it be
fig_ver = clus_count+1; % how many plots tall should it be
fspac = 0.005; % the spacing around the plots, on all sides
fpadd = 0.005; % the spacing around the plots, on all sides, this takes more space than fspac though
fmarg = 0.05; % the margins around the plots, at the edge of the figure
fsiz = 5; % the fontsize for different texts
flw = 1; % the line width for different plots

%% add an annotation to the figure with some important info
ann_str = sprintf('Tetrode: %d, Spikes: %d, Time: %d, Analysed: %s',tet,spik_count,duration,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');     

w_width = 50;
% plot original clusters
for cc = 1:clus_count
    clu = clusters(cc); % clu = the current cluster
    uci = ['r' sdata.rat_num '_' sdata.date '_t' num2str(tet) '_c' num2str(clu)]; % string for use in structure array
    cspt = sdata.(uci).spike_times;

    subaxis(fig_ver,fig_hor,cc+1,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
    [tisi,sisi] = spikeINTERVALS(cspt,w_width);
    bar(tisi,sisi,1,'k')
    xlim([-w_width/2,w_width/2])
    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
    title(sprintf('C%d',clu),'FontSize',14)

    subaxis(fig_ver,fig_hor,cc+(clus_count*cc)+1,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
    bar(tisi,sisi,1,'k')
    xlim([-w_width/2,w_width/2])
    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);
    ylabel(sprintf('C%d',clu),'FontSize',14) % label y axis
end % for cc = 1:clus_count

% plot cross-correlations
pairs = [nchoosek(clusters,2); fliplr(nchoosek(clusters,2))]; % every possible pair of clusters, including auto-correlations
for pp = 1:length(pairs(:,1))
    pnow = pairs(pp,:);
    uci1 = ['r' sdata.rat_num '_' sdata.date '_t' num2str(tet) '_c' num2str(pnow(1))]; % string for use in structure array
    uci2 = ['r' sdata.rat_num '_' sdata.date '_t' num2str(tet) '_c' num2str(pnow(2))]; % string for use in structure array

    spt1 = sdata.(uci1).spike_times;
    spt2 = sdata.(uci2).spike_times;
    [tisi,sisi] = spikeINTERVALS(spt1,w_width,spt2);
    sindx = sub2ind([clus_count+1,clus_count+1],pnow(1)+1,pnow(2)+1);

    subaxis(fig_ver,fig_hor,sindx,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
    bar(tisi,sisi,1,'b');
    xlim([-w_width/2,w_width/2]);
    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
    ax = gca;
    ax.XTick = []; % change Xtick locations to these values
    ax.YTick = []; % change Xtick locations to these values
    axis off
end % for pp = 1:length(pairs(:,1))

%% Save the figure   
figfile = [pwd 'klustest\' sdata.combined_name '\figures\'];
[~,~,~] = mkdir(figfile);
id = ['E' num2str(tet) '_cross-correlograms'];
print(fig_corr,'-dpng','-r300',[figfile id '.png']);
close(fig_corr);



































