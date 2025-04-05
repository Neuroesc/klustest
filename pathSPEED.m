 function pov = pathSPEED(pos,ppm,srate,wsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION  instantaneous speed of a trajectory
% calculate the instantaneous speed of a trajectory, given the points along that trajectory
% this will function for n-dimensional data. The pixel ratio of the recording system and the
% sampling rate can also be specified for more precise calculation
%
% USAGE:
%         pov = pathSPEED(pos,ppm,srate);
%
% INPUT:
%         pos - position data [x,y,z,d4,d5...] in pixels (or cm and set ppm to 100, or mm and set ppm to 1000)
%         ppm - pixels per metre, a single value or a vector same size as pos(:,1), default is 300
%         srate - sampling rate of the system in Hz (default is 50)
%
% OUTPUT:
%    pov - (cm/s) the velocity of the trajectory in pos, first samples will be NaN as they cannot be estimated
%
% EXAMPLES:
%     start_x = 0; 
%     start_y = 0;
%     turns = normrnd(0,30,1,1000); 
%     distances = 0.01*exprnd(5,1,1000);
%     x = start_x + cumsum(distances.*sind(cumsum(turns))); 
%     y = start_y + cumsum(distances.*cosd(cumsum(turns))); 
%     v = pathSPEED([x(:) y(:)],1000,50);
%
% See also: KLUSTEST POSTPROCESSDACQDATA

% HISTORY:
% version 1.0.0, Release 26/06/2018: Initial release
% version 1.0.0, Release 10/08/2018: Fixed bug where first sample was set to NaN instead of last, changed output to cm/s
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
if ~exist('ppm','var') || isempty(ppm)
    ppm = 300;
end
if ~exist('srate','var') || isempty(srate)
    srate = 50;
end
if ~exist('wsize','var') || isempty(wsize)
    wsize = 3; % (time samples) number specifying the window over which to estimate speed, i.e. 3 = over 3 position samples
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
p2 = circshift(pos,-1,1); % shift data one time point backward
ds = sqrt(sum((pos-p2).^2,2)); % distance between these two, gives distance travelled from previous sample at every sample
ds(end) = NaN; % last sample cannot be computed

dist_sum = movsum(ds,wsize); % moving sum of distance over window of size wsize

tsize = wsize .* (1/srate); % the length (s) of each window
pov = (dist_sum ./ tsize) ./ ppm; % speed in metres per second = distance in pixels / time in seconds / pixels per metre
pov = pov(:) .* 100; % speed in cm/s





























