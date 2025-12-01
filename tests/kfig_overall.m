function kfig_overall(sdata,uci,pname,fig_vis,save_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an sdata structure, a cell identifier and a part name and generates a figure with summary 3D plots
%   kfig_overall(sdata,uci,pname,fig_vis,save_fig)
%
%%%%%%%% Inputs
%   sdata = sdata structure
%   uci = unique cell identifier
%   pname = part name
%   fig_vis = (optional), 'on' to show figures, 'off' to hide figures, default is 'off'
%   save_fig = (optional), 1 to save figures in .fig files, 0 otherwise, default is 0
%
%%%%%%%% Comments
%   30/03/17 created to contain all of this figure code
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create figure
fig_overall = figure('visible',fig_vis,'Position',[50,60,1000,900]);
set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white
colormap(jet(256)); % to make sure the colormap is not the horrible default one
fig_hor = 3; % how many plots wide should it be
fig_ver = 3; % how many plots tall should it be
fspac = 0.01; % the spacing around the plots, on all sides
fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
fmarg = 0.03; % the margins around the plots, at the edge of the figure
fsiz = 10; % the fontsize for different texts
flw = 2; % the line width for different plots

% add an annotation to the figure with some important info
ann_str = sprintf('Cell: %s, Part: %s, Analysed: %s',uci,pname,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spikes and position plot                        
subaxis(fig_ver,fig_hor,1,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
% get data
pox = sdata.(pname).pox;
poy = sdata.(pname).poy;
poz = sdata.(pname).poz;
spx = sdata.(uci).(pname).spx;
spy = sdata.(uci).(pname).spy;
spz = sdata.(uci).(pname).spz;
lx = sdata.(pname).lx;
ly = sdata.(pname).ly;
lz = sdata.(pname).lz;
part_duration = sdata.(pname).duration;
pfrate = sdata.(uci).(pname).frate;

% plot data
plot3(pox,poy,poz,'k')
hold on
plot3(spx,spy,spz,'r.','MarkerSize',20)
h = line(lx,ly,lz);
set(h,'Color',[0.5 0.5 0.5 0.5],'LineWidth',1,'LineStyle','-')
axis off
view(3)
daspect([1 1 1])
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
axis vis3d
camproj perspective
rotate3d on

% add annotation with spatial information etc
anns = sprintf('Spk:%d\nt:%ds\nFr:%.2f',numel(spx),round(part_duration),pfrate);                  
annotation('textbox',[0.01, 0.9, 1, 0],'string',anns,'FontSize',fsiz,'LineStyle','none','interpreter','none'); 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike distribution     
subaxis(fig_ver,fig_hor,4,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
% add outlines of a cube
xpatch = [0 1 1 0 0 0; 1 1 0 0 1 1; 1 1 0 0 1 1; 0 1 1 0 0 0]*(max(pox)-min(pox))+min(pox);
ypatch = [0 0 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 0 1 1 1 1]*(max(poy)-min(poy))+min(poy);
zpatch = [0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 1 1 0 1; 1 1 1 1 0 1]*(max(poz)-min(poz))+min(poz);
for vv = 1:6
    h = patch(xpatch(:,vv),ypatch(:,vv),zpatch(:,vv),'w');
    set(h,'edgecolor','k','FaceColor','none')
end % for vv = 1:6     
hold on 

% position data
plot3(ones(size(pox)).*(max(pox)-min(pox))+min(pox),poy,poz,'Color',[0.5 0.5 0.5 0.5])
plot3(pox,poy,zeros(size(pox))+min(poz),'k','Color',[0.5 0.5 0.5 0.5])
plot3(pox,ones(size(poy)).*(max(poy)-min(poy))+min(poy),poz,'k','Color',[0.5 0.5 0.5 0.5])
% spike data
plot3(ones(size(spx)).*(max(pox)-min(pox))+min(pox),spy,spz,'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',5)
plot3(spx,spy,zeros(size(spz))+min(poz),'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',5)
plot3(spx,ones(size(spy)).*(max(poy)-min(poy))+min(poy),spz,'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',5)

axis off
view(3)
daspect([1 1 1])
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
axis vis3d
camproj perspective
rotate3d on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ratemap               
subaxis(fig_ver,fig_hor,2,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
% get data
ratemap = sdata.(uci).(pname).ratemapA;
skaggs = sdata.(uci).(part_now).spatial_measures.spatial_information;
spars = sdata.(uci).(part_now).spatial_measures.sparsity;
MMF = sdata.(uci).(part_now).spatial_measures.mean_method_focus;
lx2 = sdata.(pname).lxrm;
ly2 = sdata.(pname).lyrm;
lz2 = sdata.(pname).lzrm;

% plot data
ratemap(ratemap < 0.2*nanmax(ratemap(:))) = NaN;
vol3d('cdata',ratemap,'texture','3D');
alphamap('rampup');
% alphamap(0.5 .* alphamap);
colormap('jet');
caxis([0 nanmax(ratemap(:))]);
hold on
h = line(lx2,ly2,lz2);
set(h,'Color',[0.5 0.5 0.5 0.5],'LineWidth',1,'LineStyle','-');
axis on
view(3)
daspect([1 1 1])
xlabel('X');
ylabel('Y');
zlabel('Z');
axis vis3d
camproj perspective
rotate3d on
%title(sprintf('spars=%.1f,skaggs=%.2f,entropy=%.2f',spars,spati,entrpy),'FontSize',fsiz)

% add annotation with spatial information etc
anns = sprintf('Mx:%.1f\nMn:%.1f\nSI:%.2f\nSy:%.2f\nMMF:%.2f',nanmax(ratemap(:)),nanmean(ratemap(:)),skaggs,spars,MMF);                  
annotation('textbox',[0.31, 0.9, 1, 0],'string',anns,'FontSize',fsiz,'LineStyle','none','interpreter','none');                          
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% place field convexhulls                        
subaxis(fig_ver,fig_hor,3,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
% get data
dwellmap = sdata.(pname).dwellmapA;
fields = sdata.(uci).(pname).nfields;
fdata = sdata.(uci).(pname).field_data;
            
% plot overall position data outline
dwellmap2 = dwellmap;
dwellmap2(dwellmap < 1) = NaN; % remove voxels with a coverage less than 1 second (the minimum when generating ratemap)
[p2,p1,p3] = ind2sub(size(dwellmap2),find(~isnan(dwellmap2))); % find the subscript indices of these elements
K = convhull(p1,p2,p3); % find the convexhull of these points
tri = trisurf(K,p1,p2,p3); % plot it
set(tri,'FaceColor','k','EdgeAlpha',0.1)
alpha(.1)
hold on

% plot place field convexhulls
annotation('textbox',[0.93, 0.98, 1, 0],'string',sprintf('fields:%d',fields),'FontSize',fsiz,'LineStyle','none','interpreter','none'); 
ypos = 0.96;                        

for ff = 1:fields % for every detected place field
    % plot outline of this place field
    ps = fdata(ff).voxel_index;
    K = convhull(ps(:,2),ps(:,1),ps(:,3));
    tri = trisurf(K,ps(:,2),ps(:,1),ps(:,3));
    set(tri,'FaceColor','r');
    alpha(.1);  

    % add annotation to the plot showing the skaggs SI content of each plane
    anns = sprintf('%.2fbins\necc:%.2f\nelg:%.2f',fdata(ff).volume_voxels,fdata(ff).eccentricity,fdata(ff).elongation1);    
    annotation('textbox',[0.93, ypos, 1, 0],'string',anns,'FontSize',fsiz,'LineStyle','none','interpreter','none');                              
    ypos = ypos - 0.07;
end % for ff = 1:fields % for every detected place field                                                
axis on
view(3)
daspect([1 1 1])
xlabel('X');
ylabel('Y');
zlabel('Z');
axis vis3d
camproj perspective
rotate3d on            
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% place field convexhull projections                    
subaxis(fig_ver,fig_hor,6,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
[p2,p1,p3] = ind2sub(size(dwellmap2),find(~isnan(dwellmap2))); % find the subscript indices of these elements
K = convhull(p1,p2,p3); % find the convexhull of these points
n = length(p1);
dx = ones(n,1).*size(ratemap,1);
dy = ones(n,1).*size(ratemap,2);
dz = zeros(n,1);
px = patch('Faces',K,'Vertices',[dx p1 p3]);
py = patch('Faces',K,'Vertices',[p2 dy p3]);
pz = patch('Faces',K,'Vertices',[p2 p1 dz]);
set([px,py,pz],'FaceColor','k','FaceAlpha',.5,'EdgeColor','none');

for ff = 1:fields % for every detected place field
    ps = fdata(ff).voxel_index;
    K = convhull(ps(:,2),ps(:,1),ps(:,3));
    n = length(ps(:,1));
    dx = ones(n,1).*size(ratemap,1)-1;
    dy = ones(n,1).*size(ratemap,2)-1;
    dz = zeros(n,1)+1;
    px = patch('Faces',K,'Vertices',[dx ps(:,1) ps(:,3)]);
    py = patch('Faces',K,'Vertices',[ps(:,2) dy ps(:,3)]);
    pz = patch('Faces',K,'Vertices',[ps(:,2) ps(:,1) dz]);    
    set([px,py,pz],'FaceColor','r','FaceAlpha',.5);
end % for ff = 1:fields % for every detected place field 

axis on
view(3)
daspect([1 1 1])
xlabel('X');
ylabel('Y');
zlabel('Z');
axis vis3d
camproj perspective
rotate3d on            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
% skaggs compression matrices
subaxis(fig_ver,fig_hor,5,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
% get the data
mmap1 = sdata.(uci).(pname).spatinfo_maps{1};           
mmap2 = sdata.(uci).(pname).spatinfo_maps{2};            
mmap3 = sdata.(uci).(pname).spatinfo_maps{3};            
skag1 = sdata.(uci).(pname).spatinfo_dims(1);
skag2 = sdata.(uci).(pname).spatinfo_dims(2);
skag3 = sdata.(uci).(pname).spatinfo_dims(3);

% get the sizes of the autocorrelation maps
mmap1b = mmap2;
mmap2b = fliplr(mmap1);
mmap3b = rot90(mmap3,3);
[s1a,s2a] = size(mmap1b); % back right wall
[s1b,s2b] = size(mmap2b); % back left wall                
[s1c,s2c] = size(mmap3b); % floor       

% add outlines of a cube
xpatch = ([0 1 1 0 0 0; 1 1 0 0 1 1; 1 1 0 0 1 1; 0 1 1 0 0 0]*(s1c-1))+1;
ypatch = ([0 0 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 0 1 1 1 1]*(s2c-1))+1;
zpatch = ([0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 1 1 0 1; 1 1 1 1 0 1]*(s1a-1))+1;
for vv = 1:6
    h = patch(xpatch(:,vv),ypatch(:,vv),zpatch(:,vv),'w');
    set(h,'edgecolor','k','FaceColor','none')
end % for vv = 1:6     
hold on                        

% plot the compressed maps onto the sides of the cube                        
surf([s1c s1c; s2a s2a],[s2c 1; s2c 1],[s1a s1a; 1 1],'CData',mmap1b,'FaceColor','texturemap'); % back right wall
surf([1 s1c; 1 s1c],[s2c s2c; s2c s2c],[s1b s1b; 1 1],'CData',mmap2b,'FaceColor','texturemap'); % back left wall    
surf([1 s1c; 1 s1c],[s2c s2c; 1 1],[1 1; 1 1],'CData',mmap3b,'FaceColor','texturemap'); % floor   
axis on
view(3)
daspect([1 1 1])
xlabel('X');
ylabel('Y');
zlabel('Z');
axis vis3d
camproj perspective
rotate3d on

% add annotation to the plot showing the skaggs SI content of each plane
anns = sprintf('SI\nYZ:%.2f\nXZ:%.2f\nXY:%.2f',skag1,skag2,skag3);                  
annotation('textbox',[0.31, 0.6, 1, 0],'string',anns,'FontSize',fsiz,'LineStyle','none','interpreter','none');                     
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
% autocorrelation compression matrices            
subaxis(fig_ver,fig_hor,8,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
% get the data
auto1 = sdata.(uci).(pname).gridscore_maps{1};            
auto2 = sdata.(uci).(pname).gridscore_maps{2};            
auto3 = sdata.(uci).(pname).gridscore_maps{3};            
gscor1 = sdata.(uci).(pname).gscore_comp(1);
gscor2 = sdata.(uci).(pname).gscore_comp(2);
gscor3 = sdata.(uci).(pname).gscore_comp(3);            
            
% get the sizes of the autocorrelation maps
auto1b = auto2;
auto2b = fliplr(auto1);
auto3b = rot90(auto3,3);                        
[s1a,s2a] = size(auto1b); % back right wall
[s1b,s2b] = size(auto2b); % back left wall                
[s1c,s2c] = size(auto3b); % floor    

% add outlines of a cube
xpatch = ([0 1 1 0 0 0; 1 1 0 0 1 1; 1 1 0 0 1 1; 0 1 1 0 0 0]*(s1c-1))+1;
ypatch = ([0 0 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 0 1 1 1 1]*(s2c-1))+1;
zpatch = ([0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 1 1 0 1; 1 1 1 1 0 1]*(s1a-1))+1;
for vv = 1:6
    h = patch(xpatch(:,vv),ypatch(:,vv),zpatch(:,vv),'w');
    set(h,'edgecolor','k','FaceColor','none')
end % for vv = 1:6          
hold on

% plot the autocorrelation maps onto the sides of the cube
surf([s1c s1c; s2a s2a],[s2c 1; s2c 1],[s1a s1a; 1 1],'CData',auto1b,'FaceColor','texturemap'); % back right wall
surf([1 s1c; 1 s1c],[s2c s2c; s2c s2c],[s1b s1b; 1 1],'CData',auto2b,'FaceColor','texturemap'); % back left wall    
surf([1 s1c; 1 s1c],[s2c s2c; 1 1],[1 1; 1 1],'CData',auto3b,'FaceColor','texturemap'); % floor         
axis on
view(3)
daspect([1 1 1])
xlabel('X');
ylabel('Y');
zlabel('Z');
axis vis3d
camproj perspective
rotate3d on

% add annotation to the plot showing the grid score for each plane
anns = sprintf('Gscore\nYZ:%.2f\nXZ:%.2f\nXY:%.2f',gscor1,gscor2,gscor3);                  
annotation('textbox',[0.31, 0.5, 1, 0],'string',anns,'FontSize',fsiz,'LineStyle','none','interpreter','none');                            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D autocorrelation
subaxis(fig_ver,fig_hor,9,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
ratemap2 = ratemap;
ratemap2(isnan(ratemap)) = 0;
amap = ifftn(fftn(ratemap2).*conj(fftn(ratemap2)),'symmetric');
% amap(amap < 0.2*nanmax(amap(:))) = NaN;
vol3d('cdata',amap,'texture','3D');
alphamap('rampup');
axis on
view(3)
daspect([1 1 1])
xlabel('X');
ylabel('Y');
zlabel('Z');
axis vis3d
camproj perspective
rotate3d on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Waveforms
subplot('Position',[0.05,0.05,0.2,0.25])
wavtime = -200:20:780;
maxwavs = sdata.(uci).(pname).waveform_max;
[~,mval] = nanmax(maxwavs);
[hl,hp] = boundedline(wavtime,sdata.(uci).(pname).waveform_mean{mval},sdata.(uci).(pname).waveform_stdv{mval},'-k');
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
title(sprintf('Peak=%.1f,Width=%.f',sdata.(uci).(pname).waveform_max(mval),sdata.(uci).(pname).waveform_width(mval)),'FontSize',8);
grid on

yvec = [0.24, 0.175, 0.11, 0.045];
for ww = 1:4
    subplot('Position',[0.255,yvec(ww),0.0625,0.0625])
    wavtime = -200:20:780;
    
    [hl,hp] = boundedline(wavtime,sdata.(uci).(pname).waveform_mean{ww},sdata.(uci).(pname).waveform_stdv{ww},'-k');
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
end % for ww = 1:4

% Save the figure    
[~,~,~] = mkdir([pwd 'klustest3\' sdata.combined_name '\figures\']);
print(fig_overall,'-dpng','-r150',[pwd 'klustest3\' sdata.combined_name '\figures\' uci pname '.png'])
if save_fig
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')');
    savefig(fig_overall,[pwd 'klustest3\' sdata.combined_name '\figures\' uci pname '.fig'],'compact');
end % if save_fig
close(fig_overall);  






















































