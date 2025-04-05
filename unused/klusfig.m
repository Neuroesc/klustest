function klusfig(tet,clu,mtint,sdata,uci,part_now,fig_vis,save_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an sdata structure, a cell identifier and a part name and generates a figure with useful plots
%   klusfig(sdata,uci,part_now,fig_vis,save_fig)
%
%%%%%%%% Inputs
%   sdata = sdata structure
%   uci = unique cell identifier
%   pname = part name
%   fig_vis = (optional), 'on' to show figures, 'off' to hide figures, default is 'off'
%   save_fig = (optional), 1 to save figures in .fig files, 0 otherwise, default is 0
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

if ~exist('save_fig','var') || isempty(save_fig)
    save_fig = 0;
end % if ~exist('save_fig','var') || ismepty(save_fig)

pox = sdata.(part_now).pox;
poy = sdata.(part_now).poy;
pot = sdata.(part_now).pot;
spx = sdata.(uci).(part_now).spx;
spy = sdata.(uci).(part_now).spy;
ratemap = sdata.(uci).(part_now).ratemap;
dwellmap = sdata.(uci).(part_now).dwellmap;
automap = sdata.(uci).(part_now).grid_autocorrelation;
part_duration = sdata.(part_now).duration;
fieldd = sdata.(uci).(part_now).field_data;
skaggs = sdata.(uci).(part_now).spatial_measures.spatial_information;
spars = sdata.(uci).(part_now).spatial_measures.sparsity;
MMF = sdata.(uci).(part_now).spatial_measures.mean_method_focus;
grid_score = sdata.(uci).(part_now).grid_score;
grid_spacing = sdata.(uci).(part_now).grid_spacing;
field_size = sdata.(uci).(part_now).grid_field_size;
grid_orientation = sdata.(uci).(part_now).grid_orientation;
grid_ellipticity = sdata.(uci).(part_now).grid_ellipticity;

hd1 = sdata.(uci).(part_now).hd_session;
hd3 = sdata.(uci).(part_now).hd_cell;
hd_c = sdata.(uci).(part_now).hd_frate;
rayleigh = sdata.(uci).(part_now).hd_rayleigh;
mx2 = sdata.(uci).(part_now).hd_maximum;
ai = linspace(0,2*pi,sdata.config.hd_bins)'; % angles for binning

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create figure
fig_overall = figure('visible',fig_vis,'Position',[100, 100, 1200, 800]);
set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white
colormap(jet(256)); % to make sure the colormap is not the horrible default one
fig_hor = 4; % how many plots wide should it be
fig_ver = 3; % how many plots tall should it be
fspac = 0.01; % the spacing around the plots, on all sides
fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
fmarg = 0.03; % the margins around the plots, at the edge of the figure
fsiz = 8; % the fontsize for different texts
flw = 1; % the line width for different plots

% add an annotation to the figure with some important info
ann_str = sprintf('Cell: %s, Part: %s, Analysed: %s',uci,part_now,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');  

xv = [0 0.01 0.25 0.5 0.75];
yv = [0.04 0.16 0.43 0.7];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spike plot with black lines for path and red dots for spikes
axes('Position',[xv(2),yv(4),0.25,0.25])
plot(pox,poy,'k')
hold on
plot(spx,spy,'ro','MarkerFaceColor','r','MarkerSize',2)
daspect([1 1 1])
axis xy off
axis([min(pox)-15 max(pox)+15 min(poy)-15 max(poy)+15]);
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
title(sprintf('%d spikes (%.2f Hz)',numel(spx),numel(spx)/part_duration),'FontSize',6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dwell time heatmap
axes('Position',[xv(2),yv(3),0.25,0.25])
im = imagesc(dwellmap);
set(im,'alphadata',logical(dwellmap));
daspect([1 1 1])
axis xy off
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
title(sprintf('%.2fs (%.2f mins)',part_duration,part_duration/60),'FontSize',6);  
tt = text(-1,5,sprintf('Mx: %.2f (s), Mu: %.2f (s), Mn: %.2f (s)',nanmax(dwellmap(:)),nanmean(dwellmap(:)),nanmin(dwellmap(:))),'FontSize',6);
set(tt,'rotation',90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% firing rate map
axes('Position',[xv(3),yv(4),0.25,0.25])
im = imagesc(ratemap);
set(im,'alphadata',~isnan(ratemap));
title('Ratemap')
daspect([1 1 1])
axis xy off
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);   
% colorbar
title(sprintf('Skaggs: %.2f, Sparsity: %.2f, MMF: %.2f',skaggs,(spars*100),MMF),'FontSize',6);
tt = text(-1,5,sprintf('Mx: %.2f (Hz), Mu: %.2f (Hz), Mn: %.2f (Hz)',nanmax(ratemap(:)),nanmean(ratemap(:)),nanmin(ratemap(:))),'FontSize',6);
set(tt,'rotation',90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Place field map    
axpf = axes('Position',[xv(4),yv(4),0.25,0.25]);
bmap = fieldd.binary_ratemap;
im = imagesc(bmap);
set(im,'alphadata',~isnan(bmap));
title(sprintf('Fields: %d',length(fieldd.fields(:,1))),'FontSize',6);
daspect([1 1 1])
axis xy on
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
colormap(axpf,'bone')
ax = gca;
ax.XTick = [];
ax.YTick = [];
size_text = [];
for ff = 1:length(fieldd.fields(:,1))
    text(fieldd.fields(ff,2),fieldd.fields(ff,3),sprintf('%d',ff),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
    size_now = fieldd.fields(ff,6);
end % for ff = 1:length(fieldd.fields(:,1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% grid autocorrelation
axgs = axes('Position',[xv(5),yv(4),0.25,0.25]);
im = imagesc(automap);
set(im,'alphadata',~isnan(automap));
daspect([1 1 1])
axis xy off
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
title(sprintf('G: %.2f, Spacing: %.2f, Orientation: %.2f, Ellipticity: %.2f',grid_score,grid_spacing,grid_orientation,grid_ellipticity),'FontSize',6);
colormap(axgs,(hsv(256)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% head direction polar plot
axes('Position',[xv(5),yv(3),0.25,0.25])
polar(ai,hd1);
hold on
polar(ai,hd3);
% mmp = mmpolar([ai(:); ai(1)],hd1(:),'k:',[ai(:); ai(1)],hd3(:),'b-','FontSize',fsiz,'Grid','on','RGridVisible','off','RTickVisible','off','TTickDelta',20,'RTickLabelVisible','on','TTickLabelVisible','on');
% set(mmp,'LineWidth',1.5)
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  
title(sprintf('r: %.2f, max: %.2f',rayleigh,mx2),'FontSize',6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spikes vs time plot
axst = axes('Position',[xv(2)+0.04,yv(1),0.93,0.07]);
bstvals = 0:sdata.config.time_bins:sdata.session_duration; % vector of time points at which we should calculate spike probability
[bspikes,~] = histc(sdata.(uci).spike_times,bstvals);
bspikes = imresize(bspikes',[10,ceil(sdata.session_duration)],'Method','bilinear');
imagesc(bspikes)
colormap(axst,flipud(bone(256)))
ax = gca;
ax.YTick = [];
ylabel('Spikes')
hold on
axis xy
ax = axis;
xlabel('Time (s)')
line([min(pot) max(pot)],[0.5 0.5],'Color','r','LineWidth',4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phase preference plot     
phprobs = sdata.(uci).(part_now).spike_phase_ksdensity(:,1);
xiph = sdata.(uci).(part_now).spike_phase_ksdensity(:,2);
swav = sdata.(uci).(part_now).spike_phase_ideal;
phmu = sdata.(uci).(part_now).phase_mean;
phmud = rad2deg(phmu);
yi = sdata.(uci).(part_now).spike_phase_binned;

axes('Position',[xv(3),yv(2),0.2,0.25])
yyaxis left
ai = -pi:0.1:3*pi;
yi = [yi' yi'];
ai = ai(:);
yi = yi(:);
bar(ai,yi);
hold on
plot(ai,[swav(:); swav(:)],'b:');
% ylabel('Frequency') % label y axis
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  
ax = gca;
ax.YTick = [];

yyaxis right
plot(xiph,phprobs,'k')
title('Theta phase')
% ylabel('Probability') % label y axis
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  

xlabel('Theta Phase') % label x axis
set(gca,'Xlim',[0,2*pi]);
v = axis;
set(gca,'Ylim',[0,v(4)]);
ax = gca;
ax.XTick = linspace(0,2*pi,7); % change Xtick locations to these values
ax.XTickLabel = {'0','60','120','180','240','300','360'}; % change Xtick labels to these values
title(sprintf('Mean: %.2f (%.f)',phmu,phmud),'FontSize',6);
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Waveform plot
axes('Position',[0.05,yv(2),0.15,0.25])
wavtime = -200:20:780;
maxwavs = sdata.(uci).(part_now).waveform_max;
[~,mval] = nanmax(maxwavs);
[hl,hp] = boundedline(wavtime,sdata.(uci).(part_now).waveform_mean{mval},sdata.(uci).(part_now).waveform_stdv{mval},'-k');
set(hl,'Color','r') % line color
set(hl,'LineStyle','-') % line style
set(hl,'LineWidth',1) % line width
set(hp,'FaceColor','b') % color of area
set(hp,'FaceAlpha',1) % transparency of area
ax = gca;
ax.XLim = [-200 780];
box on
xlabel('Time (ms)')
ylabel('Amplitude (uV)')
title(sprintf('Peak: %.1f, Width: %.f',sdata.(uci).(part_now).waveform_max(mval),sdata.(uci).(part_now).waveform_width(mval)),'FontSize',6);
grid on

wpos = [0.195,0.287+0.06,0.0625,0.0625; 0.195,0.223+0.06,0.0625,0.0625; 0.195,0.158+0.06,0.0625,0.0625; 0.195,0.092+0.06,0.0625,0.0625];
for ww = 1:4
    axes('Position',wpos(ww,:))
    wavtime = -200:20:780;
    [hl,hp] = boundedline(wavtime,sdata.(uci).(part_now).waveform_mean{ww},sdata.(uci).(part_now).waveform_stdv{ww},'-k');
    set(hl,'Color','r') % line color
    set(hl,'LineStyle','-') % line style
    set(hl,'LineWidth',1) % line width
    set(hp,'FaceColor','b') % color of area
    set(hp,'FaceAlpha',1) % transparency of area
    ax = gca;
    ax.XLim = [-200 780];    
    ax.XTick = [];
    ax.YTick = [];
    box on   
    axis square
end % for ww = 1:4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spike autocorrelation - theta plot           
axes('Position',[xv(4)-0.01,yv(2),0.23,0.25])
tms1 = sdata.(uci).(part_now).theta_train(:,2);
corrdata1 = sdata.(uci).(part_now).theta_train(:,1);
thetaR = sdata.(uci).(part_now).theta_r;
thetaIndx = sdata.(uci).(part_now).theta_index;
thetaPowr = sdata.(uci).(part_now).theta_power;
thetaLin = sdata.(uci).(part_now).theta_line;

bar(tms1,corrdata1,0.9,'k');
v1 = axis;
axis([tms1(1) tms1(end) v1(3) v1(4)]);
xlabel('Time lag (ms)') % label x axis
ylabel('Probability') % label y axis
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
% set(gcf,'visible','off'); % trying to stop the figure popping up
hold on
plot(tms1,thetaLin,'r','linewidth',2);
title(sprintf('r: %.2f, i: %.2f, p: %.2f',thetaR,thetaIndx,thetaPowr),'FontSize',6);
% axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spike autocorrelation - refractory period plot   
axes('Position',[xv(5),yv(2),0.23,0.25])
tms2 = sdata.(uci).(part_now).refractory_period(:,2);
corrdata2 = sdata.(uci).(part_now).refractory_period(:,1);
RPV = sdata.(uci).(part_now).refractory_violations; % add data to structure
cont_bounds = sdata.(uci).(part_now).refractory_contamination;
censored_estimate = sdata.(uci).(part_now).censored_estimate;
tau_r = sdata.tau_r;
tau_c = sdata.tau_c;

bar(tms2,corrdata2,0.9,'k');
v1 = axis;
axis([tms2(1) tms2(end) v1(3) v1(4)]);
xlabel('Time lag (ms)') % label x axis
% ylabel('Probability') % label y axis
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);                 
% set(gcf,'visible','off'); % trying to stop the figure popping up
title(sprintf('RPV: %d, RPVc: %.2f, Censored: %.2f',RPV,cont_bounds(1),censored_estimate),'FontSize',6);
hold on
plot([-tau_r; -tau_r],[0 v(4)],'r','LineWidth',1)
plot([tau_r; tau_r],[0 v(4)],'r','LineWidth',1) 
% axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mahalanobis distance, cluster quality plot
axes('Position',[xv(3)+0.01,yv(3)+0.05,0.17,0.2])
nd = sdata.(uci).cluster_nds;
cd = sdata.(uci).cluster_cds;
isod = sdata.(uci).cluster_isolation;
lratio = sdata.(uci).cluster_lratio;
        
if ~isempty(cd) && ~isempty(nd)
    % get maximum x value
    max_d = ceil(max([max(nd);max(cd)]));
    if (isnan(max_d)) % to cope with missing channel data
        max_d = 1;
    end % if (isnan(max_d))

    xi = linspace(0,max_d,10000); % vector of values where we want to estimate ksdensity
    [vals1,~,~] = ksdensity(cd,xi,'Support',[0 max_d],'Kernel','epanechnikov','Function','pdf','Bandwidth',0.05);
    [vals2,~,~] = ksdensity(nd,xi,'Support',[0 max_d],'Kernel','epanechnikov','Function','pdf','Bandwidth',0.05);

    a1 = area(xi,vals1);
    set(a1,'FaceColor','b')
    alpha(.5)
    hold on
    a2 = area(xi,vals2);
    set(a2,'FaceColor','k')
    alpha(.5)

    set(gca,'XScale','log');
    set(gca,'Xlim',[0,max_d]);
    xlabel('Log(Distance)') % label x axis
    ylabel('Probability') % label y axis
    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
    v = axis;
    text(0.2,v(4)*0.95,sprintf('IsoD: %.2f, Lratio: %.2f',isod,lratio),'FontSize',6)
end % if ~isempty(cd) && ~isempty(nd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster space plot  
axes('Position',[xv(4)-0.01,yv(3)+0.05,0.17+0.06,0.2])
cname = sdata.combined_name;
fdata = sdata.fdata;
nfets = sdata.nfets;

clusters_now = unique(mtint.tetrode(tet).cut);
clus_count = numel(clusters_now);
clus_cut = mtint.tetrode(tet).cut;
load([pwd '\kwiktint\' cname '.kk'],'-mat','kkfet','kkset');
fetNames = kkfet.names;
fetStr = kkfet.string;
fs = find(fetStr == 1);

[~,dindx] = sort(maxwavs,2,'descend'); % find the order of these values
mch1 = dindx(1); % the channel with the maximum amplitude
mch2 = dindx(2); % the channel with the second maximum amplitude                    
plot_features = 1; % plot the first feature used (should be first principle component)

d1 = fdata(:,((mch1-1)*nfets)+plot_features(1)); % get this feature data for this channel
d2 = fdata(:,((mch2-1)*nfets)+plot_features(1)); % get this feature data for this channel
for c_plot = 1:clus_count % for every cluster on this tetrode
    cnow = clusters_now(c_plot);
    cindx = find(clus_cut == cnow);
    hold on
    colplot = [0.5 0.5 0.5 0.5];
    plot(d2(cindx),d1(cindx),'.','MarkerSize',3,'color',colplot); % plot the clusters
end % for cc = 1:clus_count
cindx = find(clus_cut == clu);
plot(d2(cindx),d1(cindx),'.','MarkerSize',3,'color','r'); % re-plot the current cluster in red to make sure it is on top and stands out

axis xy on
xlabel([fetNames{fs(plot_features(1))} ' (ch ' num2str(mch2) ')']) % label x axis
ylabel([fetNames{fs(plot_features(1))} ' (ch ' num2str(mch1) ')']) % label y axis
v = axis;
text(0.2,v(4)*0.95,sprintf('%d clusters',clus_count),'FontSize',fsiz)
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the figure    
figfile = [pwd 'klustest\' sdata.combined_name '\figures\'];
[~,~,~] = mkdir(figfile);
print(fig_overall,'-dpng','-r150',[figfile uci '_' part_now '.png'])
if save_fig
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')');
    savefig(fig_overall,[figfile uci '_' part_now '.fig'],'compact');
end % if save_fig
close(fig_overall);                                          

















































