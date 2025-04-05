function figCLUST(sdata,uci,part_now,fig_vis,save_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an sdata structure, a cell identifier and a part name and generates a figure with useful plots
%   figCLUST(tet,clu,sdata,uci,part_now,fig_vis,save_fig)
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

% determine cluster and tetrode
sst = strsplit(uci,'_');
tet = str2double(sst{3}(regexp(sst{3},'\d')));
tetstr = ['t' num2str(tet)];
clu = str2double(sst{4}(regexp(sst{4},'\d')));

% spatial data
pox = sdata.(part_now).pox;
poy = sdata.(part_now).poy;
pot = sdata.(part_now).pot;
spx = sdata.(uci).(part_now).spx;
spy = sdata.(uci).(part_now).spy;
spt = sdata.(uci).(part_now).spt;

ratemap = double(sdata.(uci).(part_now).ratemap);
dwellmap = double(sdata.(uci).(part_now).dwellmap);
automap = double(sdata.(uci).(part_now).grid_autocorrelation);
part_duration = sdata.(part_now).duration;
fieldd = sdata.(uci).(part_now).field_data;
skaggs = sdata.(uci).(part_now).spatial_measures.spatial_information;
spars = sdata.(uci).(part_now).spatial_measures.sparsity;
MMF = sdata.(uci).(part_now).spatial_measures.mean_method_focus;
MI = sdata.(uci).(part_now).spatial_measures.mutual_info;
% MIp = sdata.(uci).(part_now).spatial_measures.mutual_info_p;
overD = sdata.(uci).(part_now).over_dispersion;
grid_score = sdata.(uci).(part_now).grid_score;
grid_wavelength = sdata.(uci).(part_now).grid_metrics.wavelength;
grid_orientation = sdata.(uci).(part_now).grid_metrics.orientation;

% HD data
hd_type = sdata.config.hd_type;
if strcmp(hd_type,'histogram')
    hd1 = sdata.(uci).(part_now).hd_session;
    hd3 = sdata.(uci).(part_now).hd_cell;
    hd_c = sdata.(uci).(part_now).hd_frate;
    rayleigh = sdata.(uci).(part_now).hd_rayleigh;
    mx2 = sdata.(uci).(part_now).hd_maximum;
    mn2 = sdata.(uci).(part_now).hd_mean; % add data to structure                
    sd2 = sdata.(uci).(part_now).hd_stdev; % add data to structure                
else
    hd1 = sdata.(uci).(part_now).hd_density_session;
    hd3 = sdata.(uci).(part_now).hd_density_cell;
    hd_c = sdata.(uci).(part_now).hd_density_frate;
    rayleigh = sdata.(uci).(part_now).hd_density_rayleigh;
    mx2 = sdata.(uci).(part_now).hd_density_maximum;
    mn2 = sdata.(uci).(part_now).hd_density_mean; % add data to structure                
    sd2 = sdata.(uci).(part_now).hd_density_stdev; % add data to structure                
end
ai = linspace(0,2*pi,sdata.config.hd_bins)'; % angles for binning    

% cluster space data
nd = single(sdata.(tetstr).fetdata.noise_dists{clu});
cd = single(sdata.(tetstr).fetdata.clust_dists{clu});
isod = sdata.(tetstr).fetdata.isolation_distances(clu);
lratio = sdata.(tetstr).fetdata.lratios(clu);
fetdata = sdata.(tetstr).fetdata;
fdata = fetdata.fetdata;
nfets = fetdata.nfeatures;
clusters_now = unique(fetdata.cluster_labels);
clus_count = numel(clusters_now);
clus_cut = fetdata.cluster_labels;
tspikes = fetdata.nspikes;

% speed data
sscore = sdata.(uci).(part_now).speed_score;
sslope = sdata.(uci).(part_now).speed_slope;
sintpt = sdata.(uci).(part_now).speed_intercept;
scurve = sdata.(uci).(part_now).speed_curve;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create figure
fig_overall = figure('visible',fig_vis,'Units','pixels','Position',[100, 100, 1600, 800]);
set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white
colormap(jet(256)); % to make sure the colormap is not the horrible default one
fsiz = 8; % the fontsize for different texts
fsiz2 = 6;
flw = 1; % the line width for different plots

% add an annotation to the figure with some important info
ann_str = sprintf('Cell: %s, Part: %s, Analysed: %s',uci,part_now,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');  

xv = [0 20 240 460 680 900 1120 1340];
yv = [24 130 350 560 780 1000 1220];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spike plot with black lines for path and red dots for spikes
tvals = sdata.part_config.(part_now).times;
axes('Units','pixels','Position',[xv(3),yv(4),200,200])
hold on
for tt = 1:length(tvals(:,1))
    tindx = pot > tvals(tt,1) & pot < tvals(tt,2);
    plot(pox(tindx),poy(tindx),'k')
end
% plot spikes after so they are all on top
for tt = 1:length(tvals(:,1))
    sindx = spt > tvals(tt,1) & spt < tvals(tt,2);
    plot(spx(sindx),spy(sindx),'ro','MarkerFaceColor','r','MarkerSize',2)
end
daspect([1 1 1])
axis xy off
axis([min(pox)-15 max(pox)+15 min(poy)-15 max(poy)+15]);
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
title(sprintf('%d spikes (%.2f Hz)',numel(spx),numel(spx)/part_duration),'FontSize',6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dwell time heatmap
axdt = axes('Units','pixels','Position',[xv(2),yv(4),200,200]);
set(fig_overall,'CurrentAxes',axdt);
im = imagesc(dwellmap);
set(im,'alphadata',logical(dwellmap));
daspect([1 1 1])
axis xy off
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
title(sprintf('%.2fs (%.2f mins)',part_duration,part_duration/60),'FontSize',6);  
tt = text(-1,5,sprintf('%c: %.2f (s), %c: %.2f (s), %c: %.2f (s)',char(708),nanmax(dwellmap(:)),char(181),nanmean(dwellmap(:)),char(709),nanmin(dwellmap(:))),'FontSize',6);
set(tt,'rotation',90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% firing rate map
axrt = axes('Units','pixels','Position',[xv(4),yv(4),200,200]);
set(fig_overall,'CurrentAxes',axrt);
im = imagesc(ratemap);
set(im,'alphadata',~isnan(ratemap));
title('Ratemap')
daspect([1 1 1])
axis xy off
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);   
% colorbar
title(sprintf('SI: %.2f, Sp: %.2f, MI: %.2f, MMF: %.2f',skaggs,spars*100,MI,MMF),'FontSize',6);
tt = text(-1,5,sprintf('%c: %.2f (Hz), %c: %.2f (Hz), %c: %.2f (Hz)',char(708),nanmax(ratemap(:)),char(181),nanmean(ratemap(:)),char(709),nanmin(ratemap(:))),'FontSize',6);
set(tt,'rotation',90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Place field map    
axpf = axes('Units','pixels','Position',[xv(5),yv(4),200,200]);
set(fig_overall,'CurrentAxes',axpf);
bmap = fieldd.binary_ratemap;
im = imagesc(bmap);
set(im,'alphadata',~isnan(bmap));
daspect([1 1 1])
axis xy on
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);  
title(sprintf('%d fields, %carea: %.1fpx, %cSNR: %.1f, OD: %.1f',length(fieldd.fields(:,1)),char(181),nanmean(fieldd.fields(:,6)),char(181),nanmean(fieldd.fields(:,9)),overD),'FontSize',6);
colormap(axpf,'bone')
ax = gca;
ax.XTick = [];
ax.YTick = [];
for ff = 1:length(fieldd.fields(:,1))
    text(fieldd.fields(ff,2),fieldd.fields(ff,3),sprintf('%d',ff),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% grid autocorrelation
axgs = axes('Units','pixels','Position',[xv(6),yv(4),200,200]);
set(fig_overall,'CurrentAxes',axgs);
% cmask = sdata.(uci).(part_now).grid_metrics.central_mask;
% near_peaks = sdata.(uci).(part_now).grid_metrics.near_peaks;
% spokes = sdata.(uci).(part_now).grid_metrics.spokes;
% im = imagesc(automap);
% set(im,'alphadata',cmask);    
% hold on
% if ~isnan(spokes)
%     plot(near_peaks(:,1),near_peaks(:,2),'Marker','x','MarkerSize',10,'Color',[0.5 0.5 0.5],'LineStyle','none');
%     hold on
%     dd = drawLine(spokes);
%     set(dd,'Color',[0.5 0.5 0.5]);
% end
% daspect([1 1 1])
% axis xy off
% set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
% title(sprintf('G: %.2f, Wavelength: %.2f, Orientation: %.2f',grid_score,grid_wavelength,grid_orientation),'FontSize',6);
% colormap(axgs,(jet(256)))
msk = sdata.(uci).(part_now).grid_metrics.ring_mask;
imc = imagesc(automap);
if ~isnan(msk)
    msk(~msk) = 0.2;
    set(imc,'alphadata',msk);    
end
title(sprintf('G: %.2f, Wavelength: %.2f, Orientation: %.2f',grid_score,grid_wavelength,grid_orientation),'FontSize',6);
caxis([0 nanmax(automap(:))])
daspect([1 1 1]);
axis xy off
set(gca,'LineWidth',flw,'layer','top','FontSize',6);             
colormap(axgs,(jet(256)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% firing rate by speed
axvs = axes('Units','pixels','Position',[xv(8),yv(4)+30,260,170]);
set(fig_overall,'CurrentAxes',axvs);
plot(scurve(:,1),scurve(:,2),'k');
hold on
rr = refline(sslope,sintpt);
set(rr,'Color','r');
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
title(sprintf('Sscore: %.2f, slope: %.2f, intercept: %.2f',sscore,sslope,sintpt),'FontSize',6);
ax = gca;
ax.XLim = [0 50];
ax.YLim(1) = 0;
xlabel('Speed (cm/s)')
ylabel('Firing Rate (Hz)')
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spikes vs time plot
axst = axes('Units','pixels','Position',[xv(2)+10,yv(1),1550,60]);
set(fig_overall,'CurrentAxes',axst);
bstvals = 0:1:sdata.session_duration; % vector of 1s time points at which we should calculate spike probability
spiketime = sdata.(uci).spike_times;
[bspikes,~] = histc(spiketime,bstvals);
im = imagesc(bspikes');
colormap(axst,flipud(bone(256)))
ax = gca;
ax.YTick = [];
ylabel('Spikes')
hold on
axis xy
ax = axis;
xlabel('Time (s)')
tvals = sdata.part_config.(part_now).times;
for tt = 1:length(tvals(:,1))
    line([tvals(tt,1) tvals(tt,2)],[0.5 0.5],'Color','r','LineWidth',4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phase preference plot     
axphe = axes('Units','pixels','Position',[xv(2),yv(3)+30,200,170]);
set(fig_overall,'CurrentAxes',axphe);

phprobs = sdata.(uci).(part_now).spike_phase_ksdensity(:,1);
swav = sdata.(uci).(part_now).spike_phase_ideal;
phmu = sdata.(uci).(part_now).phase_mean;
phmud = rad2deg(phmu);
phmx = sdata.(uci).(part_now).phase_max;
phmxd = rad2deg(phmx);
phr = sdata.(uci).(part_now).phase_r;
php = sdata.(uci).(part_now).phase_p;
yi = sdata.(uci).(part_now).spike_phase_binned;

yyaxis left
aip = rad2deg(-pi:0.1:3*pi);
bar(aip,yi,1,'k');
hold on
plot(aip,swav,'b:');
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  
ax = gca;
ax.YTick = [];
ax.YColor = 'k';

yyaxis right
plot(aip,phprobs,'r')
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % Put axis lines on top, change their width and fontsize, make sure the text is black  
ax = gca;
ax.YColor = 'k';
xlabel('Theta Phase') % label x axis
set(gca,'Xlim',[min(aip) max(aip)]);
v = axis;
set(gca,'Ylim',[0,v(4)]);
astr = sprintf('r: %.1f, p: %.2f, %c: %.2f (%.f%c), %c: %.2f (%.f%c)',phr,php,char(708),phmx,phmxd,char(176),char(956),phmu,phmud,char(176));
annotation(fig_overall,'textbox','Units','pixels','Position',[xv(2)+10,yv(3)+194,600,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none');  
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Waveform plot
axes('Units','pixels','Position',[xv(1)+50,yv(2),200,200])
wavtime = -200:20:780;
maxwavs = sdata.(uci).(part_now).waveform_max;
[ax_high,mval] = nanmax(maxwavs);
ax_high = ax_high + max(max(sdata.(uci).(part_now).waveform_stdv{mval}));
ax_low = sdata.(uci).(part_now).waveform_min(mval) - max(max(sdata.(uci).(part_now).waveform_stdv{mval}));
[hl,hp] = boundedline(wavtime,sdata.(uci).(part_now).waveform_mean{mval},sdata.(uci).(part_now).waveform_stdv{mval},'-k');
set(hl,'Color','r') % line color
set(hl,'LineStyle','-') % line style
set(hl,'LineWidth',1) % line width
set(hp,'FaceColor','b') % color of area
set(hp,'FaceAlpha',1) % transparency of area
ax = gca;
ax.XLim = [-200 780];
ax.YLim = [ax_low ax_high];
box on
xlabel('Time (ms)')
ylabel(sprintf('Amplitude (%cV)',char(956)))
astr = sprintf('Peak: %.1f%cV\n Width: %.fms\n SNR: %.2f',sdata.(uci).(part_now).waveform_max(mval),char(956),sdata.(uci).(part_now).waveform_width(mval),sdata.(uci).(part_now).channel_snr(mval));
annotation(fig_overall,'textbox','Units','pixels','Position',[xv(1)+150,yv(2)+194,100,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none','HorizontalAlignment','right');      

grid on

wpos = [xv(1)+260,yv(2)+150,50,50; xv(1)+260,yv(2)+95,50,50; xv(1)+260,yv(2)+40,50,50; xv(1)+260,yv(2)-15,50,50];
for ww = 1:4
    axes('Units','pixels','Position',wpos(ww,:))
    wavtime = -200:20:780;
    
    if ~all(sdata.(uci).(part_now).waveform_mean{ww}) % if the waveform is just zeros
        ax = gca;
        ax.XLim = [-200 780]; 
        ax.YLim = [ax_low ax_high];    
        ax.XTick = [];
        ax.YTick = [];
        box on   
        axis square        
        text(-100,0,'Grounded','FontSize',7)
        continue
    end % if ~all(sdata.(uci).(part_now).waveform_mean{ww})
    
    [hl,hp] = boundedline(wavtime,sdata.(uci).(part_now).waveform_mean{ww},sdata.(uci).(part_now).waveform_stdv{ww},'-k');
    
    % we want the maximum waveform to be different
    if ww == mval
        set(hl,'Color','r') % line color
        set(hl,'LineStyle','-') % line style
        set(hl,'LineWidth',1) % line width
        set(hp,'FaceColor','b') % color of area
        set(hp,'FaceAlpha',1) % transparency of area        
    else
        set(hl,'Color','w') % line color
        set(hl,'LineStyle','-') % line style
        set(hl,'LineWidth',1) % line width
        set(hp,'FaceColor','k') % color of area
        set(hp,'FaceAlpha',1) % transparency of area
    end % if ww == mval
    
    ax = gca;
    ax.XLim = [-200 780]; 
    ax.YLim = [ax_low ax_high];    
    ax.XTick = [];
    ax.YTick = [];
    box on   
    axis square
end % for ww = 1:4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spike autocorrelation - theta plot           
axes('Units','pixels','Position',[xv(3)+20,yv(3)+30,300,170])
theta_data = sdata.(uci).(part_now).theta_train_long;
corrdata = theta_data(:,1);
tms = theta_data(:,2);
corrfilt = theta_data(:,3);
theta_index = sdata.(uci).(part_now).theta_index;
theta_ratio = sdata.(uci).(part_now).theta_ratio;
theta_skip = sdata.(uci).(part_now).theta_skip;
burst_index = sdata.(uci).(part_now).burst_index;
burst_length_mean = sdata.(uci).(part_now).burst_length_mean;

yyaxis right
plot(tms,corrfilt,'Color','r','linewidth',1);
ax = gca;
ax.YTick = [];
ax.YColor = 'k';

yyaxis left
bar(tms,corrdata,1,'k');
v1 = axis;
ax = gca;
ax.YColor = 'k';
axis([tms(1) tms(end) v1(3) v1(4)]);
xlabel('Time lag (ms)') % label x axis
ylabel('Probability') % label y axis
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
astr = sprintf('%c index: %.2f, %c ratio: %.2f, %c skip: %.2f',char(952),theta_index,char(952),theta_ratio,char(952),theta_skip);
annotation(fig_overall,'textbox','Units','pixels','Position',[xv(3)+20,yv(3)+194,600,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interspike intervals
axes('Units','pixels','Position',[xv(6)+40,yv(3)+30,300,170])

yyaxis left
adist = sdata.(uci).(part_now).isi_data.adist; % actual isi histogram   
if ~all(isnan(adist))
    bar(adist(:,1),adist(:,2),1,'k');
    hold on
    bar(-adist(:,1),adist(:,2),1,'k');
    ax = gca;
    ax.YTick = [];

    yyaxis right
    fdist = sdata.(uci).(part_now).isi_data.fdist; % fitted isi histogram   
    fx = fdist(:,1);
    plot(fdist(:,1),fdist(:,2),'r','LineStyle','-');
    hold on
    plot(-fdist(:,1),fdist(:,2),'r','LineStyle','-');

    ax = gca;
    ax.YTick = [];
    hmax = sdata.(uci).(part_now).isi_data.half_max;
    ps = sdata.(uci).(part_now).isi_data.hwidth_ps;
    line([fx(ps(2)) fx(ps(3))],[hmax hmax],'Color','r','LineWidth',1.5)
    line([fx(ps(1)) fx(ps(1))],ax.YLim,'Color','r','LineWidth',1)

    box on
    xlabel('Interspike interval (ms)') % label x axis
    set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); 
    astr = sprintf('fwhm: %.2f\nsd: %.2f',sdata.(uci).(part_now).isi_data.fwhmx,sdata.(uci).(part_now).isi_data.stdev);
    annotation(fig_overall,'textbox','Units','pixels','Position',[xv(6)+40,yv(3)+194,600,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none');  
    astr = sprintf('RPVs: %d\nRPVp: %.1f\n fpRate: %.2f',sdata.(uci).(part_now).rpv_total,sdata.(uci).(part_now).rpv_proportion,sdata.(uci).(part_now).rpv_false_positive1);
    annotation(fig_overall,'textbox','Units','pixels','Position',[xv(6)+240,yv(3)+194,100,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none','HorizontalAlignment','right');  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D phase precession
% axes('Units','pixels','Position',[xv(8)-50,yv(3)+30,270,170])
% pdata = sdata.(uci).(part_now).phase_precession;
% slope = pdata.regression_info.slope;
% phase = pdata.regression_info.phase;
% rmse = pdata.regression_info.fit_rmse;
% pmap = pdata.regression_info.img;
% res = pdata.regression_info.img_res;
% img_slope = pdata.regression_info.img_slope;
% img_phase = pdata.regression_info.img_phase;
% img_phase_offset = pdata.regression_info.img_phase_offset;
% 
% pmap = imgaussfilt(pmap,5);
% imagesc([-res/2,res/2],[-res/2,res/2],pmap);
% axis xy;   
% 
% hold on
% rline = refline(img_slope,img_phase);
% rline.Color = 'k';    
% rline.LineWidth = 2;    
% rline = refline(img_slope,img_phase+img_phase_offset);
% rline.Color = 'k';   
% rline.LineWidth = 2;    
% rline = refline(img_slope,img_phase+img_phase_offset);
% rline.Color = 'k';   
% rline.LineWidth = 2;    
% rline = refline(img_slope,img_phase-img_phase_offset);
% rline.Color = 'k'; 
% rline.LineWidth = 2;  
% axis([-res/2 res/2 -res/2 res/2])
% 
% ax.XTick = [0 res*0.25 res*0.5 res*0.75 res];
% ax.XTickLabel = {'-1','-0.5','0','0.5','1'};
% ax.YTick = [0 res*0.25 res*0.5 res*0.75 res];
% ax.YTickLabel = {'0','180','360','540','720'};
% xlabel('pdcd');
% ylabel(sprintf('%c phase',char(952)));
% set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); 
% 
% astr = sprintf('slope: %.2f\nitrcpt: %.2f\nrmse: %.2f',slope,phase,rmse);
% annotation(fig_overall,'textbox','Units','pixels','Position',[xv(8)-50,yv(3)+194,600,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none','Color','w');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spike autocorrelation - refractory period plot   
axes('Units','pixels','Position',[xv(5)-80,yv(3)+30,300,170])
tms2 = sdata.(uci).(part_now).refractory_period(:,2);
corrdata2 = sdata.(uci).(part_now).refractory_period(:,1);
tau_r = sdata.tau_r;

bar(tms2,corrdata2,0.9,'k');
v1 = axis;
axis([tms2(1) tms2(end) v1(3) v1(4)]);
xlabel('Time lag (ms)') % label x axis
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);     

astr = sprintf('burst index: %.2f',burst_index);
annotation(fig_overall,'textbox','Units','pixels','Position',[xv(5)-80,yv(3)+194,600,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none'); 
astr = sprintf('burst %clength: %.2f',char(956),burst_length_mean);
annotation(fig_overall,'textbox','Units','pixels','Position',[xv(5)+140,yv(3)+194,600,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none'); 

hold on
plot([-tau_r; -tau_r],[0 v(4)],'r','LineWidth',1)
plot([tau_r; tau_r],[0 v(4)],'r','LineWidth',1) 
plot([-6; -6],[0 v(4)],'r:','LineWidth',0.5)
plot([6; 6],[0 v(4)],'r:','LineWidth',0.5) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mahalanobis distance, cluster quality plot
axes('Units','pixels','Position',[xv(4)-100,yv(2),240,200]);
if ~isempty(cd) && ~isempty(nd)
    % get maximum x value
    max_d = ceil(max([max(nd);max(cd)]))+1;
    if (isnan(max_d)) % to cope with missing channel data
        max_d = 1;
    end 

    xi = linspace(0,max_d,10000); % vector of values where we want to estimate ksdensity
    [vals1,~,~] = ksdensity(cd,xi,'Support',[-1 max_d],'Kernel','epanechnikov','Function','pdf','Bandwidth',0.05);
    [vals2,~,~] = ksdensity(nd,xi,'Support',[-1 max_d],'Kernel','epanechnikov','Function','pdf','Bandwidth',0.05);

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
    astr = sprintf('IsoD: %.2f\nLratio: %.2f',isod,lratio);
    annotation(fig_overall,'textbox','Units','pixels','Position',[xv(4)-100,yv(2)+194,600,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none');      
    
else
    text(0.1,0.5,'Not computed');
end 
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster space plot  
axes('Units','pixels','Position',[xv(5)-50,yv(2),250,200])

% find the 2 channels with the maximum waveforms
[~,dindx] = sort(maxwavs,2,'descend'); % find the order of these values
mch1 = dindx(1); % the channel with the maximum amplitude
mch2 = dindx(2); % the channel with the second maximum amplitude                    
start_cols = 1:nfets:nfets*4; % the starting column for each channel

% get the feature data for the two channels
d1 = fdata(:,start_cols(mch1)); % get this feature data for this channel
d2 = fdata(:,start_cols(mch2)); % get this feature data for this channel

% plot the noise in grey and then the cell in red
colplot = [0.5 0.5 0.5 0.5];
plot(d1,d2,'.','MarkerSize',3,'color',colplot); % plot the clusters
hold on
cindx = find(clus_cut == clu);
plot(d1(cindx),d2(cindx),'.','MarkerSize',3,'color','r'); % re-plot the current cluster in red to make sure it is on top and stands out

axis xy on
astr = ['Ch' num2str(mch1) ' vs ' 'Ch' num2str(mch2)];
annotation(fig_overall,'textbox','Units','pixels','Position',[xv(5)-50,yv(2)+194,600,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none');  
astr = sprintf('%d clusters, %d spikes (%.1f%% of total %d)',clus_count,numel(spx),numel(spx)/tspikes*100,tspikes);
annotation(fig_overall,'textbox','Units','pixels','Position',[xv(5)-50,yv(2)+6,600,10],'string',astr,'FontSize',fsiz2,'LineStyle','none','interpreter','none');  
set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% head direction polar plot
axhd = axes('Units','pixels','Position',[xv(7),yv(4),200,200]);
set(fig_overall,'CurrentAxes',axhd);
mmp = mmpolar(ai,hd1,'k',ai,hd3,'b','FontSize',fsiz,'Grid','on','RGridVisible','off','RTickVisible','off','TTickDelta',20,'RTickLabelVisible','on','TTickLabelVisible','on');
set(mmp(1),'LineWidth',0.5);
set(mmp(2),'LineWidth',0.5);
p1 = patch(get(mmp(1),'XData'),get(mmp(1),'YData'),'k','FaceAlpha',0.2);
p2 = patch(get(mmp(2),'XData'),get(mmp(2),'YData'),'b','FaceAlpha',0.5);
title(sprintf('r: %.2f, %c: %.2f, %c: %.2f, %s: %.2f',rayleigh,char(708),mx2,char(956),mn2,char(963),sd2),'FontSize',6);
hl = legend([p1,p2],{'Session','Cell'},'Position',[xv(3)+0.12,yv(3),0.05,0.05]);
legend boxoff
set(hl,'FontSize',6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the figure    
figfile = [pwd '\klustest\' sdata.combined_name '\figures\'];
[~,~,~] = mkdir(figfile);
print(fig_overall,'-dpng','-r150',[figfile uci '_' part_now '.png'])
if save_fig
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')');
    savefig(fig_overall,[figfile uci '_' part_now '.fig'],'compact');
end % if save_fig
% close(fig_overall);                                          

















































