function [hd_smoo] = estimateHD(pox,poy,srate,twindow,smoo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A short function to estimate an animal's heading direction based on position data, hd_smoo is the estimated HD in degrees
% pox = x position coordinates
% poy = y position coordinates
% srate = sampling rate of the position data in ms
% twindow (optional) = the length of window to estimate HD over in seconds (default = 1s)
% smoo (optional) = smoothing factor in data points for estimated HD (default = 5)
%
% 25/02/16 Created
% 01/03/16 Fixed, added twindow, fixed angles so they deviate between [0,360], added defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('twindow','var') || isempty(twindow) 
    twindow = srate;
else
    twindow = twindow * srate;
end % if ~exist('smoo','var') || isempty(smoo) 
if ~exist('smoo','var') || isempty(smoo) 
    smoo = 5;
end % if ~exist('smoo','var') || isempty(smoo)

%% Estimate heading at every time point
hd_smoo = NaN(size(pox));
for p = 1:length(pox)-twindow
    x1 = pox(p);
    x2 = pox(p+twindow);
    y1 = poy(p);
    y2 = poy(p+twindow);

    ang_deg = atan2(y2-y1,x2-x1)*180/pi;
    hd_smoo(p) = ang_deg; 
end % for p = 1:length(pox)-1

%% Figure for checking
fig_hd = figure('visible','off');
rowa1 = 2; % Number of rows for subplots
colu1 = 2; % Number of columns for subplots
sa_spac1 = 0.03; % Spacing between plots
sa_padd1 = 0.03;
sa_marg1 = 0.01;
fsize1 = 4;
axlwid = 1;

subaxis(rowa1,colu1,1,'Spacing',sa_spac1,'Padding',sa_padd1,'Margin',sa_marg1);
plot(pox,poy,'k')
hold on
plot(pox(1),poy(1),'ro')
daspect([1 1 1])

subaxis(rowa1,colu1,2,'Spacing',sa_spac1,'Padding',sa_padd1,'Margin',sa_marg1);
rose(deg2rad(hd_smoo),30)

subaxis(rowa1,colu1,3,'Spacing',sa_spac1,'Padding',sa_padd1,'Margin',sa_marg1);
scatter(pox,poy,15,hd_smoo)
colormap([jet(128); flipud(jet(128))])
daspect([1 1 1])
colorbar

subaxis(rowa1,colu1,4,'Spacing',sa_spac1,'Padding',sa_padd1,'Margin',sa_marg1);
binz = linspace(-180,180,60);
binned = histc(hd_smoo,binz);
binned = binned/max(binned);
bar(binz,binned)
axis([-180 180 0 1])

antext = ['Session estimated head direction | Time window: ' num2str(twindow) ' | Smoothing factor: ' num2str(smoo)];
annotation('textbox', [0.05, 1.0, 1.0, 0], 'string',antext,'FontSize',5);                    % [x y w h]                                                         
[~,~,~] = mkdir('Figures');
saveas(fig_hd,['Figures\Estimated_HD_summary.png'],'png');
close(fig_hd);


        
        
        
        
        
        
        
        
        
        
        
        
        
        





