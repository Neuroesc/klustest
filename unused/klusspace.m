function klusspace(tet,mtint,sdata,fig_vis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function generates a figure showing the cluster space of a tetrode
%   klusspace(tet,mtint,sdata,fig_vis)
%
%%%%%%%% Inputs
%   tet = the tetrode to look at
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

cname = sdata.combined_name;
clusters_now = unique(mtint.tetrode(tet).cut);
clusters_now(clusters_now == 0) = [];
clus_count = numel(clusters_now);
spik_count = sdata.spike_count{tet};
duration = sdata.session_duration;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create figure
fig_clu = figure('visible',fig_vis,'Position',[100, 100, 1024, 800]);

set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white
colormap(jet(256)); % to make sure the colormap is not the horrible default one
fig_hor = 4; % how many plots wide should it be
fig_ver = 2; % how many plots tall should it be
fspac = 0.03; % the spacing around the plots, on all sides
fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
fmarg = 0.03; % the margins around the plots, at the edge of the figure
fsiz = 5; % the fontsize for different texts
flw = 1; % the line width for different plots

%% add an annotation to the figure with some important info
ann_str = sprintf('Session: %s, Tetrode: %d, Spikes: %d, Time: %d, Analysed: %s',cname,tet,spik_count,duration,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');     

%% cluster space
[~,~,~,~,fdata,nfets,~] = clusterQUALITY(cname,tet); % get the session feature data for this tetrode
clusters = unique(mtint.tetrode(tet).cut);
clus_count = numel(clusters);
clus_cut = mtint.tetrode(tet).cut;
load([pwd '\kwiktint\' cname '.kk'],'-mat','kkfet','kkset');
fetNames = kkfet.names;
fetStr = kkfet.string;
fs = find(fetStr == 1);

% I was going to add this as a user argument, but the features are not arranged in a convenient way (i.e. this actually means plot the first and second included - which may be
% any two features. It doesn't mean plot features 1 and 2, which would be PC1 and PC2. In the .fet file the features are just concatenated with no way to tell which is which
% I might come back to make this option more flexible - the used features are in fetNames and fetStr, the features themselves are in fdata
plot_features = [1]; 

annotation('textbox',[0.55, 1, 1, 0],'string',fetNames{fs(plot_features(1))},'FontSize',15,'LineStyle','none','interpreter','none');     
epairs = nchoosek(1:4,2); % every possible combination of channel pair  
for pp = 1:length(epairs)           
    subaxis(fig_ver,fig_hor,pp,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
    c1 = epairs(pp,1); % first channel to plot (x axis)
    c2 = epairs(pp,2); % second channel to plot (y axis)
    d1 = fdata(:,((c1-1)*nfets)+plot_features(1));
    d2 = fdata(:,((c2-1)*nfets)+plot_features(1));

    colmap = jet(clus_count);
    linfo = cell(1,clus_count);
    for cc_plot = 1:clus_count
        cnow = clusters(cc_plot);
        cindx = find(clus_cut == cnow);
        hold on
        plot(d2(cindx),d1(cindx),'.','MarkerSize',2,'color',colmap(cc_plot,:));
        axis xy on square
        title([num2str(c1) ' vs ' num2str(c2)]);
        xlabel(fetNames{fs(plot_features(1))}); % label x axis
        ylabel(fetNames{fs(plot_features(1))}); % label y axis
        set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
        linfo{cc_plot} = ['Cluster ' num2str(cnow)];
    end % for cc = 1:clus_count
end % for pp = 1:length(epairs)   
warning('off','MATLAB:legend:IgnoringExtraEntries'); % the legend function will want to warn that there are too many plots (lines)
warning('off','MATLAB:handle_graphics:exceptions:SceneNode');
[legh,objh,~,~] = legend(linfo,'boxoff');
M = findobj(objh,'type','Line');
set(M,'MarkerSize',50);
set(legh,'Position',[0.5 0.03 0.14 0.44],'FontSize',14);

%% Save the figure 
figfile = [pwd 'klustest\' sdata.combined_name '\figures\'];
[~,~,~] = mkdir(figfile);
id = [cname '_E' num2str(tet) '_cluster_space'];
print(fig_clu,'-dpng','-r300',[figfile id '.png'])
close(fig_clu);




























