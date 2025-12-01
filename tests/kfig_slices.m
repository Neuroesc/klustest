function kfig_slices(sdata,uci,part_now,dtype,fig_vis,save_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an sdata structure, a cell identifier and a part name and generates a figure with sliced 3D plots
%   kfig_slices(sdata,uci,pname,fig_vis,save_fig)
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

if ~exist('dtype','var') || isempty(dtype)
    dtype = 1;
end % if ~exist('dtype','var') || ismepty(dtype)
if dtype == 1 % if we want to look at spatial information content per slice
    slices = sdata.(uci).(part_now).slices;
    xdnow = slices{1};
    ydnow = slices{2};
    zdnow = slices{3};
    values = sdata.(uci).(part_now).spatinfo_slice;
    xvnow = values{1};
    yvnow = values{2};
    zvnow = values{3};
    vname = 'SI';
    fname = 'SIslices';
elseif dtype == 2 % if we want to look at grid score per slice
    slices = sdata.(uci).(part_now).auto_slices;
    xdnow = slices{1};
    ydnow = slices{2};
    zdnow = slices{3};
    values = sdata.(uci).(part_now).gscore_slice;
    xvnow = values{1};
    yvnow = values{2};
    zvnow = values{3};
    vname = 'G';
    fname = 'Gslices';
end % if dtype == 1

% count the number of slices
ns = numel(xdnow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create figure
fig_slice = figure('visible',fig_vis,'Position',[100, 100, 1024, 800]);
set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
set(gcf,'color','w'); % makes the background colour white
colormap(jet(256)); % to make sure the colormap is not the horrible default one
fig_hor = ns; % how many plots wide should it be
fig_ver = 5; % how many plots tall should it be
fspac = 0.01; % the spacing around the plots, on all sides
fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
fmarg = 0.03; % the margins around the plots, at the edge of the figure
fsiz = 8; % the fontsize for different texts
flw = 1; % the line width for different plots

% add an annotation to the figure with some important info
ann_str = sprintf('Cell: %s, Part: %s, Analysed: %s',uci,part_now,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ratemap               
subaxis(fig_ver,fig_hor,1:ns*2,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
% get data
ratemap = sdata.(uci).(part_now).ratemapA;
lx2 = sdata.(part_now).lxrm;
ly2 = sdata.(part_now).lyrm;
lz2 = sdata.(part_now).lzrm;

% plot data
ratemap(ratemap < 0.2*nanmax(ratemap(:))) = NaN;
vol3d('cdata',ratemap,'texture','3D');
alphamap('rampup');
colormap('jet');
caxis([0 nanmax(ratemap(:))]);
hold on
h = line(lx2,ly2,lz2);
set(h,'Color',[0.5 0.5 0.5 0.5],'LineWidth',1,'LineStyle','-');
axis([min(lx2)-1 max(lx2)+1 min(ly2)-1 max(ly2)+1 min(lz2)-1 max(lz2)+1])
axis on
view(3)
daspect([1 1 1])
xlabel('X');
ylabel('Y');
zlabel('Z');
axis vis3d
camproj perspective
rotate3d on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% slice plots, for each dimension
for pp = 1:length(xdnow)
    subaxis(fig_ver,fig_hor,pp+ns*2,'Spacing',0.001,'Padding',0.001,'Margin',fmarg);                        
    mnow = rot90(xdnow{pp});
    sentropy = entropy(double(mnow));    
    %mnow(isnan(mnow)) = 0;
    im = imagesc(mnow);
    set(im,'alphadata',~isnan(mnow));
    daspect([1 1 1]);
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    axis off    
    title(sprintf('%s=%.1f,SE=%.1f',vname,xvnow(pp),sentropy),'FontSize',fsiz);
    
    if pp == 1
        ylabel('Slices along Y');
    end % if pp == 1
end % for pp = 1:length(xdnow)

for pp = 1:length(ydnow)
    subaxis(fig_ver,fig_hor,pp+ns*3,'Spacing',0.001,'Padding',0.001,'Margin',fmarg);                        
    mnow = fliplr(rot90(ydnow{pp}));
    sentropy = entropy(double(mnow));    
    %mnow(isnan(mnow)) = 0;
    im = imagesc(mnow);
    set(im,'alphadata',~isnan(mnow));
    daspect([1 1 1]);
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    axis off    
    title(sprintf('%s=%.1f,SE=%.1f',vname,yvnow(pp),sentropy),'FontSize',fsiz);  
    
    if pp == 1
        ylabel('Slices along X');
    end % if pp == 1
end % for pp = 1:length(xdnow)

for pp = 1:length(zdnow)
    subaxis(fig_ver,fig_hor,pp+ns*4,'Spacing',0.001,'Padding',0.001,'Margin',fmarg);                        
    mnow = fliplr(rot90(zdnow{pp}));
    sentropy = entropy(double(mnow));    
    %mnow(isnan(mnow)) = 0;
    im = imagesc(mnow);
    set(im,'alphadata',~isnan(mnow));
    daspect([1 1 1]);
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    axis off
    title(sprintf('%s=%.1f,SE=%.1f',vname,zvnow(pp),sentropy),'FontSize',fsiz);  
    
    if pp == 1
        ylabel('Slices along Z');
    end % if pp == 1
end % for pp = 1:length(xdnow)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the figure   
figfile = [pwd 'klustest3\' sdata.combined_name '\figures\'];
[~,~,~] = mkdir(figfile);
print(fig_slice,'-dpng','-r150',[figfile uci part_now fname '.png'])
if save_fig
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')');
    savefig(fig_slice,[figfile uci part_now fname '.fig'],'compact');
end % if save_fig
close(fig_slice);                                          












































