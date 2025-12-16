function [svals,stime,sscore,sslope,sintcpt,crve,crven,crvestd] = getSPEEDmod(ppot,ppov,pspt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getSPEEDmod  calculate speed score, slope, intercept etc for a cell
% Using the analysis from Kropff, Carmichael, Moser and Moser (2015) Speed cells in the medial entorhinal cortex
% We can generate a smoothed firing rate histogram for a cell (instantaneous firing rate) and correlate this with
% instantaneous speed (speed score). We can also extract the slope and y-intercept of this relationship
%
% USAGE:
%         [stime,sscore,sslope,sintcpt,crve,crven] = getSPEEDmod(ppot,ppov,pspt)
%
% INPUT:
%         ppot - position times
%         ppov - speed at every position time
%         pspt - time of every spike
%
% OUTPUT:
%    svals - speeds for which we have values (I would suggest clipping at 50cm/s)
%    stime - time spent at every speed
%    sscore - speed score, i.e. the correlation between firing rate and running speed
%    sslope - the slope of the LLS line of best fit between firing rate and speed
%    sintcpt - y-intercept of this line
%    crve - the actual line, x values are svals (first output)
%    crven - the line normalised against the y-intercept (I never actually had to use this) x values are svals (first output)
%
% EXAMPLES:
%
% See also: klustest corr polyfit

% HISTORY:
% version 1.0.0, Release 09/04/19 Initial release, contained here to clean up klustest
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with input variables
    narginchk(3,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
    % Taken from Kropff, Carmichael, Moser and Moser (2015) Speed cells in the medial entorhinal cortex
    % Instantaneous firing rate was obtained by dividing the whole session into 20-ms bins, coinciding 
    % with the frames of the tracking camera. A temporal histogram of spiking was then obtained, smoothed 
    % with a 250-ms-wide Gaussian filter.
    [ppotu,uidx] = unique(ppot); % sometimes dacqUSB repeats time stamps at the end of a file
    edg = ppotu(1)-.01 : .02 : ppotu(end)+.01; % 20ms time bins spanning the session
    shist = histcounts(pspt,edg); % spike histogram
    fh = fspecial('gaussian',[1 13]); % gaussian filter, they say 250ms wide, but samples are every 20ms, so I have used 13 samples or 260ms instead
    shist = imfilter(shist,fh,'replicate','same'); % smooth histogram
    shist = shist ./ 0.02; % convert spikes to firing rate

    % here we bin the speed data and calculate the mean firing rate per 2 cm/s bin
    ppov = fillmissing(ppov,'nearest'); % at least one speed sample can never be computed, so we need to replace that NaN value
    ppov = ppov(uidx); % only take the values from the unique time values
    povidx = fix(ppov ./ 2)+1; % bin the speed data in 2cm/s increments
    povidx = interp1(ppotu,povidx,movmean(edg,2,'EndPoints','discard'),'nearest');
    povidx(isnan(povidx(:))) = 1000;
    A = accumarray(povidx(:),shist(:),[],@(x) mean(x,'omitnan'));
    Aerr = accumarray(povidx(:),shist(:),[],@(x) std(x,'omitnan'));
    stime = accumarray(povidx(:),ones(size(povidx))) .* (1/50); % the total time spent moving at each speed

    % The speed score for each cell was defined as the Pearson product-moment correlation between the cell’s 
    % instantaneous firing rate and the rat’s instantaneous running speed, on a scale from ?1 to 1.
    svals = ((1:length(A)).*2)' - 1; % values of speed bins in cm/s
    fvals = A(:); % mean firing rate for each bin
    sscore = corr(svals(svals<=50),fvals(svals<=50),'type','Pearson','rows','pairwise'); % speed score is the correlation between instantaneous firing rate and speed

    % Because of the variability in baseline and slope, a simple or normalized average of speed cell activity would not properly 
    % capture the population behaviour. To obtain a better normalization method, we applied to any firing rate measure f of a 
    % speed cell expressed in Hz the linear transformation [equation in paper]
    % where A (Hz) and B (cm?1) are the y intercept and slope of the cell’s speed 
    % tuning. The 50 value is given in cm s?1. This linear transformation aims to achieve for every cell a normalized dimensionless 
    % activity of 0 when the rat is still and 1 when it runs at 50 cm s?1, allowing for proper population averaging.
    P = polyfit(svals(svals<=50),fvals(svals<=50),1);
    fvals2 = (fvals - P(2)) / (P(1) * 50);

    % accumulate data
    sslope = P(1);
    sintcpt = P(2);
    crve = fvals;
    crven = fvals2;
    crvestd = Aerr(:);


    
    
%     figure
%     plot(svals,stime)
%     
%     figure
%     edg = 0:2:100;
%     f = histcounts(ppov,edg);
%     x = histcents(edg);
%     plot(x,f)
%     
%     figure
%     plot(ppov)
%     
%     keyboard



































