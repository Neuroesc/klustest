function [datout] = processPOS(pox,poy,pot,hd,ppm,cutv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datout = structure of processed data
%
% 12/03/16 created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datout = struct;
if ~exist('ppm','var') || isempty(ppm)
    ppm = 400; % pixel ratio of camera pixels/meter
end % if isempty(ppm) || ~exist(ppm)
if ~exist('cutv','var') || isempty(cutv)
    cutv = 0.5; % max velocity in ms
end % if isempty(ppm) || ~exist(ppm)
if ~exist('hd','var') || isempty(hd)
    hd = ones(size(pox)); % hd, if it is not given
end % if isempty(ppm) || ~exist(ppm)

%% Prepare position data
nindx = [find(pox == 1023 & poy == 1023); find(pox == 0 & poy == 0); find(pox == inf & poy == inf); find(pox == -inf & poy == -inf)];
pox(nindx) = NaN;
poy(nindx) = NaN;
datout.missing_values = length(find(isnan(pox) & isnan(poy)));
datout.missing_values_percent = ((length(find(isnan(pox) & isnan(poy))))/length(pox)) * 100;

%% Speed filter position data
vec1 = [pox poy; NaN NaN];
vec2 = [NaN NaN; pox poy];

[dists] = fastDist(vec1,vec2); % distance travelled at every point from last point
times = [pot; NaN] - [NaN; pot]; % time passage for every point
speeds = (dists./ppm)./times; % instantaneous speed in ms
datout.time = times(1:end-1);
datout.distance = dists(1:end-1);
datout.speed = speeds(1:end-1);

dindx = find(speeds > cutv);
for i = 1:2:(length(dindx)-1)
    pox(dindx(i):dindx(i+1)) = NaN;
    poy(dindx(i):dindx(i+1)) = NaN;
    hd(dindx(i):dindx(i+1)) = NaN;
end % for i = 1:2:(length(dindx)-1)
datout.uspox = pox;
datout.uspoy = poy;
datout.ushd = hd;

%% Smooth position data
smdata = smoothn({pox poy},10,'robust'); % a smoothing factor of 100 is equivalent to the kalman filter, but makes less of a difference as it gets larger
pox = smdata{1};
poy = smdata{2};
hd = hd_interp(hd); % linear iterpolation of missing hd data
datout.pox = pox;
datout.poy = poy;
datout.hd = hd;

%% Estimate head direction from displacement
hd2 = NaN(size(pox));
vindx = 1:(length(pox)-1);
hd2(vindx) = mod((180/pi)*(atan2(-poy(vindx+1) + poy(vindx), pox(vindx+1) - pox(vindx))),360);
hd2(end) = hd2(end-1);
datout.hd2 = hd2;

%% Calculate continuously sampled data
itimes = 0:0.02:ceil(max(pot));
datout.interp_time = itimes';
datout.interp_pox = interp1(pot,pox,itimes)';
datout.interp_poy = interp1(pot,poy,itimes)';
datout.interp_hd = interp1(pot,hd,itimes)';
datout.interp_hd2 = interp1(pot,hd2,itimes)';
datout.interp_speed = interp1(pot,speeds(1:end-1),itimes)';
datout.interp_distance = interp1(pot,dists(1:end-1),itimes)';



