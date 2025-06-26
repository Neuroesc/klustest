function [thi,thf,thr,t500,c500,f500] = getTHETAintrinsic(pspt,t500,c500)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getTHETAintrinsic  fits a decomposing sine wave to a spike autocorrelogram
% Given some spike times or a spike autocorrelogram this function fits a decomposing sine wave to the data
% in order to extract the intrinsic theta frequency and power (index) in the spike autocorrelogram
% this method is taken from: van der Meer & Redish (2011) Theta phase precession in rat ventral striatum links place and reward information.
% This method was described for use with a 500ms autocorrelogram, which is what this function will make if necessary, but it should work with other
% time lags (untested)
%
% USAGE:
%         [thi,thf,thr,t500,c500] = getTHETAintrinsic(pspt,t500,c500)
%
% INPUT:
%         pspt - vector of spike times
%         t500 - time lag values for a spike autocorrelogram
%         c500 - probability values for a spike autocorrelogram
%
%         give the function pspt and it will generate its own 500ms autocorrelogram
%         OR give it t500 and c500 if you already have them for a speed increase
%
% OUTPUT:
%    thi - theta index (power/modulation), high values = more theta present
%    thf - frequency of the detected theta in Hz
%    thr - rsquare of fit, higher means better
%    t500 - time lag values for the spike autocorrelogram used
%    c500 - probability values for the spike autocorrelogram used
%    f500 - the result of the fit at the time values in t500
%
% See also: klustest fit spikeINTERVALS_v2

% HISTORY:
% version 1.0.0, Release 09/04/19 Initial release, to contain code for klustest
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with input variables
    narginchk(1,3);
    
% preallocate outputs    
    thi = NaN; 
    thf = NaN; 
    thr = NaN;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
    % van der Meer & Redish (2011) Theta phase precession in rat ventral striatum links place and reward information
    % To quantify the degree and frequency of theta modulation in single cells, we used the method used by Royer et al. (2010). 
    % First we computed the autocorrelogram of the cell, in 10 ms bins from ?500 to +500 ms, normalized it to the maximum value 
    % between 100 and 150 ms (corresponding to theta modulation), and clipped all values above 1. Then we fit the following function:
    % (a*(sin(2*pi*w*x+(pi/2))+1)+b) .* exp(-abs(x)/t1) + (c*exp(-(x.^2)/t2.^2))
    % where t is the autocorrelogram time lag, and a–c, ?, and ?1–2 were fit using the fminsearch optimization function in MATLAB. 
    % A measure of theta modulation strength, the “theta index,” was defined as a/b, which intuitively corresponds to the ratio of 
    % the sine fit relative to the baseline in the autocorrelogram. We elected to use this measure rather than spectral methods because 
    % of the direct estimate of the theta frequency ? it provides, which in the case of a spectral estimator needs to be obtained 
    % indirectly. For best-fitting performance, we restricted possible values for ? to [6, 12], for a and b to non-negative values, 
    % for c to [0, 0.2], and for ?2 to [0, 0.05]. Accuracy of this fitting procedure was further improved by only including spike 
    % data from when animals were active (moving their heads at >3 cm/s) as described above (see Behavioral analysis) and by 
    % only fitting cells with at least 100 spikes available. Examples of the fit output can be seen in Figures 3 and 4.
    if ~exist('t500','var') || isempty(t500) || ~exist('c500','var') || isempty(c500) % if a spike autocorrelogram was not provided, make one
        [c500,t500,e] = spike_auto(pspt,'spt2',pspt,'bin_size',10,'win_size',500,'method','correlogram','log',0);
        f500 = NaN(size(t500));
    end
    if numel(pspt)<100 || isempty(pspt) || ~any(c500)
        return
    end
    
    tval = max(c500(t500>100 & t500<150),[],'omitnan');
    if tval>0 && ~isnan(tval)
        cvecn = c500 ./ tval;
        cvecn(cvecn>=1) = NaN;
        
        if sum(cvecn(:),'omitnan')==0
            return
        end

        % Fit sine wave to data  
        tms2 = t500/1000;
        nindx = isnan(tms2) | isinf(tms2) | isnan(cvecn) | isinf(cvecn);
        [xData,yData] = prepareCurveData(tms2(~nindx),cvecn(~nindx));
        ft = fittype( '(a*(sin(2*pi*w*x+(pi/2))+1)+b) .* exp(-abs(x)/t1) + (c*exp(-(x.^2)/t2.^2))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off'; opts.Lower = [0 0 0 -Inf 0 6]; opts.StartPoint = [0.2238 0.7512 0.2551 0.5059 0.6990 0.8909]; opts.Upper = [Inf Inf 0.8 Inf 0.05 12];
        try
            [fitresult,gof] = fit(xData,yData,ft,opts);
        catch
            thi = NaN;
            thf = NaN;
            thr = NaN;
            f500 = NaN(size(t500));
            return
        end
        thi = fitresult.a ./ fitresult.b;
        thf = fitresult.w;
        thr = gof.rsquare;    
        f500 = fitresult(t500/1000).*tval;
    end







































