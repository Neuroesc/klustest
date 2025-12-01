function mtint = postprocessDACQDATA(mtint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	This is a function for post-processing mtint position data, it serves as a replacement for postprocess_DACQ_data
%   and is intended to much faster. It also takes an mtint directly, rather than subfields and uses as few subfunctions
%   as possible.
%   mtint = postprocessDACQDATA(mtint)
%
%%%%%%%% Inputs
%   mtint = mtint structure, like the one made by getDACQDATA
%
%%%%%%%% Outputs
%   mtint = mtint structure, now with additional fields:
%         mtint.pos.led_pos = led positions now filtered for jumps
%         mtint.pos.jumpyPercent = the percentage of data that was found to be jumpy (speed exceeded maximum) after smoothing pos and jumps equally
%         mtint.pos.xy_pixels = animals smoothed position, either first led or mean of both leds
%         mtint.pos.xy_cm = same as above but in cm
%         mtint.pos.dir = heading direction calculated from two LEDs or from displacement if this is not available
%         mtint.pos.dir_displacement = heading direction estimated from displacement
%         mtint.pos.speed = instantaneous running speed (calculated using the smoothed, jumpless data) in cm/s
%
%%%%%%%% Comments
%   13/04/17 modified
%   13/04/17 reduced precision of outputs to reduce mtint file size
%   19/04/17 rewritten to be faster and simpler, changed name to postprocessDACQDATA, function takes mtint directly and calls few subfunctions
%   19/04/17 changed smoothing to utilise smoothn
%   19/04/17 changed jump detection for speed and simplicity, but will remove more data as smoothing window increases
%   24/07/17 added compatibility with variable pixel ratio, using getDACQDATAHEADERS
%   10/08/18 updated for klustest update, added pathSPEED
%   04/04/19 moved vectoral pixel ratio to here
%
%   © Roddy Grieves: rmgrieves@gmail.com (not sure who wrote the original postprocess_pos_data though)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial variables
jumpspeed = 1; % max speed in m/s
mtint.pos.led_pos = single(mtint.pos.led_pos); % make sure this field is in the correct precision or else NaNs become zeros
n_leds = size(mtint.pos.led_pos,2); % the number of LEDs used in recording
srate = mtint.pos.header(1).sample_rate_num;
sm_wind = round(srate/2); % window length for smoothing data and removing jumps

% pixel ratio
pratio = [];
for pr = 1:numel([mtint.pos.header.pixels_per_metre])
    pratio = [pratio; ones(mtint.pos.header(pr).num_pos_samples,1) .* mtint.pos.header(pr).pixels_per_metre];
end
mtint.pos.pixels_per_metre = uint16(pratio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process data
mtint.pos.raw_pos = mtint.pos.led_pos;

%% For 2 spot tracking, check for instances of swapping (often happens if one LED bright, one less bright).
if n_leds == 2 && ~isempty(mtint.pos.led_pix) % only check if we have 2 LEDs and the LED size (led_pix)
    swap_list = led_swap_filter(double(mtint.pos.led_pos),double(mtint.pos.led_pix));
    dum = mtint.pos.led_pos(swap_list,1,:);
    mtint.pos.led_pos(swap_list,1,:) = mtint.pos.led_pos(swap_list,2,:);
    mtint.pos.led_pos(swap_list,2,:) = dum;
    dum = mtint.pos.led_pix(swap_list, 1);
    mtint.pos.led_pix(swap_list,1) = mtint.pos.led_pix(swap_list,2);
    mtint.pos.led_pix(swap_list,2) = dum;
end

%% Speed filter points for those moving impossibly fast, replace with NaN
for ll = 1:n_leds
    % get data
    pox = double(mtint.pos.led_pos(:,ll,1));
    poy = double(mtint.pos.led_pos(:,ll,2));
    pot = (1:length(pox)) .* (1/srate);

    nindx = pox==0 & poy==0 | isnan(pox) | isnan(poy) | isinf(pox) | isinf(poy) | pox==1023 & poy==1023; % find all the stupid position data values
    pox(nindx) = NaN;
    poy(nindx) = NaN;
    
    % calculate instantaneous speed, this will have few zero values because unsmoothed tracking jumps around all the time
    pox1 = [pox(:); NaN];
    poy1 = [poy(:); NaN];
    pox2 = [NaN; pox(:)];
    poy2 = [NaN; poy(:)];
    pot1 = [pot(:); NaN];
    pot2 = [NaN; pot(:)];

    ds = sqrt((pox1-pox2).^2 + (poy1-poy2).^2); % distance between each position point and the next one
    dsa = movmean(ds,sm_wind,'omitnan','Endpoints','shrink'); % get a moving average
    ts = mean([pot1 pot2],2); % time values of the average points (they fall between each position data point)
    ds = interp1(ts(2:end-1),dsa(2:end-1),pot); % interpolate to find distance values at actual timestamps
    vs = (ds' ./ pratio) .* srate; % calculate speed in metres/second

    % find jumps in the data
    jumps = vs > jumpspeed; % find points where the speed exceeds maximum
    jumps = logical(movmean(jumps,sm_wind,'omitnan','Endpoints','shrink')); % smooth this logical vector to match the smoothing in position data (to try and catch both sides of every jump)
    jpercent = sum(jumps) ./ numel(jumps);
    mtint.pos.jumpyPercent = jpercent * 100;

    % remove them
    pox(jumps) = NaN;
    poy(jumps) = NaN;    
    
    % add data back to mtint
    mtint.pos.led_pos(:,ll,1) = single(pox);
    mtint.pos.led_pos(:,ll,2) = single(poy);  
end

%% Interpolate and smooth to replace all NaN led positions
for ll = 1:n_leds
    % smooth position data using smoothn algorithm (also interpolates NaNs)
    [datout,s,e] = smoothn({mtint.pos.led_pos(:,ll,1),mtint.pos.led_pos(:,ll,2)},100,'robust');
    mtint.pos.led_pos(:,ll,1) = datout{1};
    mtint.pos.led_pos(:,ll,2) = datout{2};       
end

%% Get heading direction
% Need to know angles (and distances) of LEDs from rat (in bearing_colour_i in .pos header)
% 2 LEDs are assumed to be at 180deg to each other with their midpoint over the
% animals head. Convention is that colbearings_set(1) is the large light (normally
% at front) and (2) is the small light. Position of lights (defined in set
% file header) is defined in degs with 0 being straight ahead of animal, values
% increase anti-clockwise.
if n_leds == 1  
    pox = mtint.pos.led_pos(:,1,1);
    poy = mtint.pos.led_pos(:,1,2);  
    mtint.pos.xy_pixels = single([pox(:) poy(:)]);
    mtint.pos.xy_cm = mtint.pos.xy_pixels ./ [pratio pratio] .* 100;
    
elseif n_leds == 2
%     colbearings_set = NaN(2,1);
%         
%     if isfield(mtint.header,'bearing_colour_1')
%         colbearings_set(1) = mtint.header(1).bearing_colour_1; % get led bearings from set file
%         colbearings_set(2) = mtint.header(1).bearing_colour_2;
%     else
%         colbearings_set(1) = mtint.header(1).lightBearing_1; % get led bearings from set file
%         colbearings_set(2) = mtint.header(1).lightBearing_2;        
%     end
    colbearings_set = zeros(2,1); % override, just assume LEDs are parallel to rat's head

    % get direction
    correction = colbearings_set(1); % to correct for light pos relative to rat subtract angle of large light   
    dir = mod((180/pi) * (atan2(-mtint.pos.led_pos(:,1,2) + mtint.pos.led_pos(:,2,2), mtint.pos.led_pos(:,1,1) - mtint.pos.led_pos(:,2,1))) - correction, 360);
    mtint.pos.dir = dir(:);

    % Get position from smoothed individual lights
    mtint.pos.xy_pixels = permute(nanmean(mtint.pos.led_pos,2),[1,3,2]);
    mtint.pos.xy_cm = mtint.pos.xy_pixels ./ pratio .* 100;
end

%% Get estimated heading from position displacement
% y from tracker increases downwards so dir is positive anticlockwise from X axis, 0<= dir <360
pox = mtint.pos.xy_pixels(:,1);
poy = mtint.pos.xy_pixels(:,2);

dird = mod((180/pi) * (atan2(-poy(2:end) + poy(1:end-1), pox(2:end) - pox(1:end-1))),360);
dird = [dird(1); dird];
dird(end) = dird(end-1);
dird = dird(:);

mtint.pos.dir_displacement = dird;
if ~isfield(mtint.pos,'dir') % if true direction is missing use displacement to always fill it
    mtint.pos.dir = dird;
end

%% Calculate instantaneous speed using smoothed data, this should contain the correct number of zero values and no values exceeding maxspeed
mtint.pos.speed = single(pathSPEED([pox(:) poy(:)],pratio,srate)); % speed in cm/s

































