function rdata = getCIRCLINregress(px,py)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes circular data (such as spike theta phase) and linear data (such as position in place field)
%   and fits 'barber' poles to it using a form of linear regression.
%   The slope of the regression is estimated in passes, with each pass the slope is refined using the sum of residuals
%   as an indicator of the best slope.
%   The phase is then estimated by testing phases around the data to again find the best one that minimises the sum
%   of residuals.
%   This analysis will not work well as the slope reaches infinity, but no analysis will anyway
%   Lastly, the circular-linear correlation is calculated using functions from the circular statistics toolbox
%   rdata = getCIRCLINregress(px,py)
%
%%%%%%%% Inputs
%   px = x (linear) predictor
%   py = p (circular) response
%
%%%%%%%% Outputs
%   rdata = structure containing results
%     rdata.slope = the slope of the best fit line
%     rdata.phase = the phase (y-intercept) of the best fit line
%     rdata.slope_deg = same as slope but in degrees
%     rdata.phase_deg = same as phase but in degrees
%     rdata.fit_rmse = the root mean squared error of the fit
%     rdata.fit_mean = the mean residual error of the fit
%     rdata.circ_corr = the circular-linear correlation [r-value,p-value]
%     rdata.img = a density map of the data
%     rdata.img_slope = same as slope but appropriate for the density map
%     rdata.img_phase = same as phase but appropriate for the density map
%     rdata.img_phase_offset = spacing between lines in the density map (equivalent to 360 degrees elsewhere)
%
%%%%%%%% Comments
%   15/11/17 created
%   15/11/17 circular regression slope working
%   15/11/17 implemented iterative approach to slope fitting
%   16/11/17 circular regression phase working (implemented differently to the paper though)
%   16/11/17 circular-linear correlation added from circular statistics toolbox
%   16/11/17 added density map of data and transformation to lines
%   16/11/17 added warning for excessive slopes
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
fit_passes = 3; % number of passes to make through data to find best slope
fit_power = 2^10; % the number of slopes to test at each pass, 2^10 is probably enough
fit_scope = 2^10; % the range of slopes to test initially, i.e. start_cent +/- fit_scope
scope_coef = 0.33; % at each pass the scope is reduced to tighten around best slope, this is the coefficient used, 0.1 is added each pass to prevent it reaching zero
start_cent = 0; % starting slope, default is flat (0)
start_phase = pi; % inbitial phase used when calculating slope, doesn't really change anything
phase_passes = 1000; % number of phase values to test when trying to find the best fit

rdata = struct;

% % to simulate data:
% ps = 3000;
% px = -1 + 2*rand(ps,1);
% py = px*0.8 + normrnd(0,1,ps,1);
% pndx = px > -1 & px < 1 & py > -pi & py < pi;
% px = px(pndx);
% py = py(pndx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate regression
%   analysis taken (as closely as possible) from Kempter  et al. (2012) Quantifying circular–linear associations: Hippocampal phase precession
%   Their phase estimation is not really understandable to me, so I have implemented something else
%   I also use an iterative approach to try and speed things up a bit
% prepare matrices
px = single(px(:));
py = single(py(:));
pxmat = ones(numel(px),fit_power) .* px;
pymat = ones(numel(py),fit_power) .* py;

%% find best slope
for pp = 1:fit_passes
%     disp(sprintf('\tpass: %d',pp));
%     disp(sprintf('\tpreferred slope: %.2f',start_cent));    
%     disp(sprintf('\tfit power: %.2f',fit_power));    
%     disp(sprintf('\tfit scope: %.2f',fit_scope));    
    
    % sort out settings for this pass
    start_slope = start_cent-fit_scope;
    end_slope = start_cent+fit_scope;
    
    % test slopes against data
    svec = single(linspace(start_slope,end_slope,fit_power)); % vector of slope values to test
    lymat = svec .* pxmat + start_phase; % matrix of y values, generated vectorally using line equation
    
    % find the 'best' one (highest vector length of residuals)
    angs = angdiff(pymat,lymat); % angle difference between actual values (py) and expected values (ly), these are our residuals
    rv = circ_r(angs,[],[],1); % mean resultant vector length R down columns of angs
    [~,rmax] = nanmax(rv); % find rmax, the rayleigh vector length with the highest value (i.e the slope with the most concentrated residuals)
    
    % make this the next start_cent
    start_cent = svec(rmax);
    best_slope = start_cent;

    % reduce fit_scope on each pass, as we get more precise we need less tests
    fit_scope = fit_scope*scope_coef+1;
    scope_coef = scope_coef+0.1;
end
% if abs(best_slope) > 9
%     disp(sprintf('\tWARNING: regression slope (%.1f) may exceed test tolerances...',best_slope));
% end

%% find phase offset
pxmat = ones(numel(px),phase_passes) .* px;
pymat = ones(numel(py),phase_passes) .* py;
smat = ones(numel(px),phase_passes) .* best_slope; % matrix of slope values

% test phases against data
% this tests a bunch of phases from the minimima to the maxima of the data, plus some margin, as long as the data cross the y-axis this should work
pvec = single(linspace(min(py)*1.9,max(py)*1.9,phase_passes)); % vector of phase values to test, I use a coefficient of 1.9 to prevent wrap around back to best phase
lymat = smat .* pxmat + pvec; % matrix of y values, generated vectorally using line equation

% find the 'best' one (lowest mean residual)
angs = angdiff(pymat,lymat); % angle difference between actual values (py) and expected values (ly), these are our residuals
mv = nansum(abs(angs),1); % sum of residual down columns of angs
[~,rmin] = nanmin(mv); % find rmin, the lowest total residual value (i.e the phase with the lowest residuals)
fit_rmse = sqrt(circ_mean(circ_std(angs(:,rmin)).^2)); % calculate fit root mean squared error (in radians, measure of how good fit is, lower is better)
fit_mean = nanmean(angs(:,rmin)); % mean residual

% get the phase value
best_phase = pvec(rmin);

% perform circular-linear correlation on data points
% taken from circular statistics toolbox: Zar (2010) Biostatistical Analysis
[r,p] = circ_corrcl(px,py);

% generate a 2D image map of relationship
py2 = [rad2deg(py)-360; rad2deg(py); rad2deg(py)+360; rad2deg(py)+720];
px2 = [px; px; px; px];
res = 128;
xq = linspace(-1,1,res);
yq = linspace(0,720,res);
[pmap,~,~,~,~] = histcounts2(py2,px2,yq,xq);
img_slope = (rad2deg(best_slope) ./ (720/res)) / (res/2);
img_phase = rad2deg(best_phase) ./ (720/res);
img_phase_offset = 360 ./ (720/res);

% accumulate data
rdata.slope = best_slope;
rdata.phase = best_phase;
rdata.slope_deg = rad2deg(best_slope);
rdata.phase_deg = rad2deg(best_phase);
rdata.fit_rmse = fit_rmse;
rdata.fit_mean = fit_mean;
rdata.circ_corr = [r,p];
rdata.img = pmap;
rdata.img_slope = img_slope;
rdata.img_phase = img_phase;
rdata.img_phase_offset = img_phase_offset;
rdata.img_res = res;



% make figures if desired
if 0
    figure
    subplot(2,2,1)
    plot(mv,pvec)
    
    
    subplot(2,2,2)
    plot(px,py,'k.')
    % axis([-1 1 0 2*pi]);
    hold on
    rline = refline(best_slope,best_phase);
    rline.Color = 'r';
    rline.LineWidth = 2;
    rline = refline(best_slope,best_phase-2*pi);
    rline.Color = 'r';
    rline.LineWidth = 2;    
    rline = refline(best_slope,best_phase+2*pi);
    rline.Color = 'r';   
    rline.LineWidth = 2;    
    rline = refline(best_slope,best_phase+4*pi);
    rline.Color = 'r';    
    rline.LineWidth = 2;    
    axis([-1 1 -pi pi])
    title(sprintf('Slope: %.1f, Intercept: %.1f, rmse: %.1f, r: %.2f, p: %.2f',best_slope,best_phase,fit_rmse,r,p));
    

    subplot(2,2,3)
    py2 = rad2deg(py);
    plot(px,py2,'k.');
    hold on
    plot(px,py2+360,'k.');
    plot(px,py2+720,'k.');
    plot(px,py2-360,'k.');
    
    best_phase2 = rad2deg(best_phase);
    best_slope2 = rad2deg(best_slope);
    rline = refline(best_slope2,best_phase2);
    rline.Color = 'r';    
    rline.LineWidth = 2;    
    rline = refline(best_slope2,best_phase2+360);
    rline.Color = 'r';   
    rline.LineWidth = 2;    
    rline = refline(best_slope2,best_phase2+720);
    rline.Color = 'r';   
    rline.LineWidth = 2;    
    rline = refline(best_slope2,best_phase2-360);
    rline.Color = 'r'; 
    rline.LineWidth = 2;    
    axis([-1 1 0 720])

    
    subplot(2,2,4)
    pmap = imgaussfilt(pmap,5);
    imagesc([-res/2,res/2],[-res/2,res/2],pmap);
    axis xy;   

    hold on
    rline = refline(img_slope,img_phase);
    rline.Color = 'k';    
    rline.LineWidth = 2;    
    rline = refline(img_slope,img_phase+img_phase_offset);
    rline.Color = 'k';   
    rline.LineWidth = 2;    
    rline = refline(img_slope,img_phase+img_phase_offset);
    rline.Color = 'k';   
    rline.LineWidth = 2;    
    rline = refline(img_slope,img_phase-img_phase_offset);
    rline.Color = 'k'; 
    rline.LineWidth = 2;  
    axis([-res/2 res/2 -res/2 res/2])

    
    
end






















