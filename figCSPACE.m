function figCSPACE(fetdata,figfile,fig_vis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function generates a figure showing the cluster space of a tetrode
%   figCSPACE(fetdata)
%
%%%%%%%% Inputs
%   tet = the tetrode to look at
%   mtint = an mtint structure
%   sdata = an sdata structure
%   fig_vis = figure visibility, set to 'on' to see figure, 'off' to hide figure
%
%%%%%%%% Comments
%   31/03/17 created to contain all of this figure code
%   05/04/17 updated figure
%   05/04/17 removed dependency on kwiktint output, now the function just takes a structure produced by clusterQUALITY
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('fig_vis','var') || isempty(fig_vis)
    fig_vis = 'off';
end % if ~exist('fig_vis','var') || ismepty(fig_vis)

maxspik = 2^18; % the maximum number of spikes to plot in cluster figures

fdata = fetdata.fetdata;
nclus = fetdata.nclusters;
nspikes = fetdata.nspikes;
nfets = fetdata.nfeatures;
isods = fetdata.isolation_distances;
lrats = fetdata.lratios;
noise_dists = fetdata.noise_dists;
clust_dists = fetdata.clust_dists;
clus = fetdata.cluster_labels;
clusters = unique(clus(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure
fig_clu = figure('visible',fig_vis,'Position',[100, 100, 1400, 800]);
set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white
colormap(jet(256)); % to make sure the colormap is not the horrible default one
fsiz = 8; % the fontsize for different texts
flw = 1; % the line width for different plots

% % add an annotation to the figure with some important info
ann_str = sprintf('Features: %d (plus time), Clusters: %d (1 noise), Spikes: %d (%d in clusters = %.1f%%), Analysed: %s',nfets,nclus,nspikes,length(fdata(clus~=0,1)),length(fdata(clus~=0,1))/nspikes*100,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none'); 

xvec = [0.04 0.04 0.04 0.28 0.52 0.76];
yvec = [0.05 0.36 0.67 0.67 0.67 0.67];
combos = nchoosek(1:nfets:nfets*4,2); % the different combinations of channels, the numbers given here actually correspond directly to columns in fdata, for first feature (probably PC1)
cols = jet(nclus);
cols(1,:) = [0.5 0.5 0.5];
totvars = zeros(6,1);
for pp = 1:6
    axes('Position',[xvec(pp),yvec(pp),0.2,0.25])
    f1 = fdata(:,combos(pp,1));
    f2 = fdata(:,combos(pp,2));    
    totvars(pp) = std(single(f1)) + std(single(f2));
    f1 = int16(f1);
    f2 = int16(f2); 
    
    % downsample data if necessary
    clus2 = clus;
    if numel(f1) > maxspik
        rindx = randi(numel(f1),[maxspik,1]);
        f1 = f1(rindx,:);
        f2 = f2(rindx,:);
        clus2 = clus(rindx,:);
    end % if numel(f1) > maxspik

    for cc = 1:nclus
        cnow = clusters(cc);
        cindx = clus2 == cnow;
        plot(f1(cindx,:),f2(cindx,:),'Marker','o','MarkerSize',3,'LineStyle','none','MarkerFaceColor',cols(cc,:),'MarkerEdgeColor','none');
        hold on
    end % for cc = 0:nclus
    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
    xlabel(sprintf('Ch %d',ceil(combos(pp,1)/4)));
    ylabel(sprintf('Ch %d',ceil(combos(pp,2)/4)));
end % for pp = 1:length(xvec)

[~,bcomb] = max(totvars);
axes('Position',[xvec(bcomb),yvec(bcomb),0.2,0.25],'Color','none')
box on
ax = gca;
ax.XTick = [];
ax.YTick = [];
set(gca,'LineWidth',flw*2); 

% downsample data if necessary
axes('Position',[0.3,0.05,0.66,0.54])
f1 = int16(fdata(:,combos(bcomb,1)));
f2 = int16(fdata(:,combos(bcomb,2))); 
if numel(f1) > maxspik
    rindx = randi(numel(f1),[maxspik,1]);
    f1 = f1(rindx,:);
    f2 = f2(rindx,:);
    clus2 = clus(rindx,:);    
end % if numel(f1) > maxspik

cnames = cell(nclus,1);
for cc = 1:nclus
    cnow = clusters(cc);    
    cindx = clus2 == cnow;
    plot(f1(cindx,:),f2(cindx,:),'Marker','o','MarkerSize',3,'LineStyle','none','MarkerFaceColor',cols(cc,:),'MarkerEdgeColor','none');
    hold on
    if cnow
        cnames{cc} = sprintf('C%d: isod=%.1f',cnow,isods(cnow));
    else
        cnames{cc} = 'Noise';        
    end % if cc
end % for cc = 0:nclus
set(gca,'LineWidth',flw*2,'layer','top','FontSize',fsiz); 
xlabel(sprintf('Ch %d',ceil(combos(bcomb,1)/4)));
ylabel(sprintf('Ch %d',ceil(combos(bcomb,2)/4)));
legend(cnames,'Location','SouthEast','Box','off')

%% Save the figure 
%print(fig_clu,'-dpng','-r300',figfile) % print can cause problems because of the large number of data points in the figure
saveas(gcf,figfile,'png');
close(fig_clu);


% old version
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Create figure
% fig_clu = figure('visible',fig_vis,'Position',[100, 100, 1024, 800]);
% 
% set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
% set(gcf,'color','w'); % makes the background colour white
% colormap(jet(256)); % to make sure the colormap is not the horrible default one
% fig_hor = 4; % how many plots wide should it be
% fig_ver = 2; % how many plots tall should it be
% fspac = 0.03; % the spacing around the plots, on all sides
% fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
% fmarg = 0.03; % the margins around the plots, at the edge of the figure
% fsiz = 5; % the fontsize for different texts
% flw = 1; % the line width for different plots
% 
% %% add an annotation to the figure with some important info
% ann_str = sprintf('Session: %s, Tetrode: %d, Spikes: %d, Time: %d, Analysed: %s',cname,tet,spik_count,duration,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
% annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');     
% 
% %% cluster space
% [~,~,~,~,fdata,nfets,~] = clusterQUALITY(cname,tet); % get the session feature data for this tetrode
% clusters = unique(mtint.tetrode(tet).cut);
% clus_count = numel(clusters);
% clus_cut = mtint.tetrode(tet).cut;
% load([pwd '\kwiktint\' cname '.kk'],'-mat','kkfet','kkset');
% fetNames = kkfet.names;
% fetStr = kkfet.string;
% fs = find(fetStr == 1);
% 
% % I was going to add this as a user argument, but the features are not arranged in a convenient way (i.e. this actually means plot the first and second included - which may be
% % any two features. It doesn't mean plot features 1 and 2, which would be PC1 and PC2. In the .fet file the features are just concatenated with no way to tell which is which
% % I might come back to make this option more flexible - the used features are in fetNames and fetStr, the features themselves are in fdata
% plot_features = [1]; 
% 
% annotation('textbox',[0.55, 1, 1, 0],'string',fetNames{fs(plot_features(1))},'FontSize',15,'LineStyle','none','interpreter','none');     
% epairs = nchoosek(1:4,2); % every possible combination of channel pair  
% for pp = 1:length(epairs)           
%     subaxis(fig_ver,fig_hor,pp,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
%     c1 = epairs(pp,1); % first channel to plot (x axis)
%     c2 = epairs(pp,2); % second channel to plot (y axis)
%     d1 = fdata(:,((c1-1)*nfets)+plot_features(1));
%     d2 = fdata(:,((c2-1)*nfets)+plot_features(1));
% 
%     colmap = jet(clus_count);
%     linfo = cell(1,clus_count);
%     for cc_plot = 1:clus_count
%         cnow = clusters(cc_plot);
%         cindx = find(clus_cut == cnow);
%         hold on
%         plot(d2(cindx),d1(cindx),'.','MarkerSize',2,'color',colmap(cc_plot,:));
%         axis xy on square
%         title([num2str(c1) ' vs ' num2str(c2)]);
%         xlabel(fetNames{fs(plot_features(1))}); % label x axis
%         ylabel(fetNames{fs(plot_features(1))}); % label y axis
%         set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
%         linfo{cc_plot} = ['Cluster ' num2str(cnow)];
%     end % for cc = 1:clus_count
% end % for pp = 1:length(epairs)   
% warning('off','MATLAB:legend:IgnoringExtraEntries'); % the legend function will want to warn that there are too many plots (lines)
% warning('off','MATLAB:handle_graphics:exceptions:SceneNode');
% [legh,objh,~,~] = legend(linfo,'boxoff');
% M = findobj(objh,'type','Line');
% set(M,'MarkerSize',50);
% set(legh,'Position',[0.5 0.03 0.14 0.44],'FontSize',14);
% 
% 






















