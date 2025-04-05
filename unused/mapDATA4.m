function [ratemap,config] = mapDATA4(pos,spk,config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function bins spike and position data into a firing rate map and can use various approaches
%   This version has been reworked so that x,y data can be easily corrected to match bins in the map. 
%   It has also been changed to accept only data in cm, not pixels.
%   [ratemap,config] = mapDATA4(pos,spk,config)
%
%%%%%%%% Inputs
%   pos = position data to make map with [x,y]
%   spk = spike data to make map with [x,y]
%   config = structure containing the settings needed when creating the map:
%         config.rmethod = (default 'KDE') the mapping approach to use, either 'nearest','gaussian','adaptive','KDE'
%         config.map_padd = (default 4) the number of bins to pad spatial maps with
%         config.bin_size = (default 2) bin size in cm for calculating the rate map (make sure data is in cm or pm/sm values are given)
%         config.map_sigma = (default 0.5) used by traditional and KDE method, sigma of gaussian to use when smoothing traditional dwell and spike maps, or used as bandwidth of KDE
%         config.min_dwell = (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
%         config.g_sigma = (default 10) only used by gaussian method - sigma of gaussian used to weight position point and spike point distances, bigger means smoother maps
%         config.min_dist = (default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
%         config.srate = (default 50) sampling rate of data in Hz, used to calculate time
%         config.limits = (default is data limits) limits of map to make, default is [min(x) max(x) min(y) max(y)]
%   pm = (optional) vector same length as pos specifying pixels per metre for each position sample, or a single integer if same pixel ratio for all
%   sm = (optional) vector same length as spk specifying pixels per metre for each spike, or a single integer if same pixel ratio for all
%
%%%%%%%% Outputs
%   ratemap = generated firing rate map
%   config = see above, settings used for map will be found here, also added are:
%         config.ratemap = the final ratemap
%         config.dwellmap = the final dwell time map, if one was created
%         config.spikemap = the final spike map if one was created
%         config.poxnew = x position data converted to match ratemap
%         config.poynew = y position data converted to match ratemap
%         config.spxnew = x spike data converted to match ratemap
%         config.spynew = y spike data converted to match ratemap
%
%%%%%%%% Comments
%   08/08/16 created as a simple mapping function
%   20/04/17 created from mapDATA because I want a function I can use with overdispersion
%   20/04/17 function now outputs converted position data (to match map)
%   20/04/17 added limit capability so spike and time maps can be made with same limits
%   27/07/17 changed to v3, this version expects data to be in cm, not pixels (meant to work with variable pixel ratio klustest)
%   01/08/17 changed to v4, overhauled, added switches for different approaches
%   02/08/17 added adaptive method
%   03/08/17 added KDE approach, comments
%   31/10/17 added better functionality for converting units of data
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('config','var')
    config = struct;
end
if ~isfield(config,'rmethod')
    config.rmethod = 'nearest';
end
if ~isfield(config,'map_padd')
    config.map_padd = 2;
end
if ~isfield(config,'bin_size')
    config.bin_size = 2;
end
if ~isfield(config,'map_sigma')
    config.map_sigma = 2;
end
if ~isfield(config,'min_dwell')
    config.min_dwell = 0.1;
end
if ~isfield(config,'g_sigma')
    config.g_sigma = 10;
end
if ~isfield(config,'min_dist')
    config.min_dist = 1;
end
if ~isfield(config,'max_dist')
    config.max_dist = config.bin_size*1.5;
end
if ~isfield(config,'srate')
    config.srate = 50;
end

% extract data
pox = double(pos(:,1));
poy = double(pos(:,2)); 
spx = double(spk(:,1));
spy = double(spk(:,2)); 

% sort out position and spike data units
if isfield(config,'pos_ppm')
    pox = pox ./ config.pos_ppm .* 100; % convert position data to cm
    poy = poy ./ config.pos_ppm .* 100; % convert position data to cm        
end
if isfield(config,'spk_ppm')
    spx = spx ./ config.spk_ppm .* 100; % convert position data to cm
    spy = spy ./ config.spk_ppm .* 100; % convert position data to cm     
end

% prepare limits of map
lx = [nanmin(pox) nanmax(pox)];
ly = [nanmin(poy) nanmax(poy)];
config.limits = [nanmin(pox) nanmax(pox) nanmin(poy) nanmax(poy)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare bins
bin_size = config.bin_size;
p = config.map_padd;
rmethod = config.rmethod;

% centre data on the origin
pox = pox - mean(lx); 
poy = poy - mean(ly);   
spx = spx - mean(lx); 
spy = spy - mean(ly);   
lx = lx - mean(lx); 
ly = ly - mean(ly);  

% vectors for bin edges (used by nearest, gaussian, adaptive methods)
xvec = 0 : bin_size : (max(abs(lx))+bin_size+bin_size*p); 
yvec = 0 : bin_size : (max(abs(ly))+bin_size+bin_size*p);    

% sort them to increase from -max to max
xvec = unique(sort([-xvec xvec],'ascend')); 
yvec = unique(sort([-yvec yvec],'ascend'));
config.xloc = xvec;
config.yloc = yvec;

% prepare bin centres (used by KDE method)
[X1,Y1] = meshgrid(xvec,yvec); 
xcen = X1(:);
ycen = Y1(:);
config.xv = xcen;
config.yv = ycen;

config.poxnew = (pox ./ bin_size) + ceil(length(xvec)/2);
config.poynew = (poy ./ bin_size) + ceil(length(yvec)/2);
config.spxnew = (spx ./ bin_size) + ceil(length(xvec)/2);
config.spynew = (spy ./ bin_size) + ceil(length(yvec)/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate map
switch rmethod
    case {1,'nearest'} % 3D histogram
        % this method is the most widely used and can be found everywhere in the literature. It is also the most basic and the fastest.
        % data are literally binned where they are found. This is done for position data and the bin counts are then multiplied by the
        % sampling rate to get time. The same is then done for spikes. The spikemap and dwellmap are then gaussian smoothed and then
        % the spikemap is divided by the dwellmap. Empty bins, or bins where the rat spent less than config.min_dwell time are empty (NaN)
        spikemap = hist3([spy,spx],{yvec,xvec});
        if config.map_sigma ~= 0
            spikemap = imgaussfilt(spikemap,config.map_sigma);
        end
        if isfield(config,'dwellmap') && ~isempty(config.dwellmap)
            dwellmap = config.dwellmap;
        else        
            dwellmap = hist3([poy,pox],{yvec,xvec}); 
            if config.map_sigma == 0
                dwellmap = dwellmap .* (1/config.srate);                
            else
                dwellmap = imgaussfilt(dwellmap,config.map_sigma) .* (1/config.srate);
            end
        end
        ratemap = spikemap ./ dwellmap;
        ratemap(dwellmap < config.min_dwell | dwellmap == 0 | isnan(dwellmap) | isinf(dwellmap)) = NaN;
        config.dwellmap = dwellmap;
        config.spikemap = spikemap;     
        config.ratemap = ratemap; 
        
    case {2,'gaussian'}
        % this is an approach taken from: Leutgeb et al. (2005) Independent Codes for Spatial and Episodic Memory in Hippocampal Neuronal Ensembles
        % for every bin we calculate the distance to every spike and every position data point. These distances are weighted with a gaussian so that
        % data close to the bin centre = 1 and data far away = 0. The sum of the spike distances are then divided by the sum of the position
        % distances to obtain firing rate. The standard deviation of the gaussian can be changed, the result being that bins become 'larger' or rather
        % bins will overlap more as more data is contained in each. This method can be slow when there are many bins.
        if ~isfield(config,'dwellmap')
            dwellmap = zeros(numel(yvec),numel(xvec));
        else
            dwellmap = config.dwellmap;
        end
        spikemap = zeros(numel(yvec),numel(xvec));

        for bb = 1:length(xcen) % for every bin   
            if ~isfield(config,'dwellmap')
                % get the distance to every position point
                rindx = pox>xcen(bb)-config.max_dist & pox<xcen(bb)+config.max_dist & poy>ycen(bb)-config.max_dist & poy<ycen(bb)+config.max_dist;            
                dp = sqrt((pox(rindx)-xcen(bb)).^2 + (poy(rindx)-ycen(bb)).^2); 
                if nanmin(dp) > config.min_dist % if this bin is too far from the position data
                    continue
                end 
                dp_norm = (exp(-0.5 * ((dp - 0)./config.g_sigma).^2) ./ (sqrt(2*pi) .* config.g_sigma)) ./ (1 ./ (sqrt(2*pi) .* config.g_sigma)); % convert the distances to gaussian weighted values
                tbin = nansum(dp_norm)*(1/config.srate); % time in voxel
                dwellmap(yvec==ycen(bb),xvec==xcen(bb)) = tbin; % accumulate data
            else
                tbin = config.dwellmap(yvec==ycen(bb),xvec==xcen(bb));
            end

            % get the distance to every spike
            if tbin < config.min_dwell % if the rat spent too little time in this voxel
                spikemap(yvec==ycen(bb),xvec==xcen(bb)) = 0;
                continue
            end
            rindx = spx>xcen(bb)-config.max_dist & spx<xcen(bb)+config.max_dist & spy>ycen(bb)-config.max_dist & spy<ycen(bb)+config.max_dist;            
            ds = sqrt((spx(rindx)-xcen(bb)).^2 + (spy(rindx)-ycen(bb)).^2); % get the distance to every spike
            ds_norm = (exp(-0.5 * ((ds - 0)./config.g_sigma).^2) ./ (sqrt(2*pi) .* config.g_sigma)) ./ (1 ./ (sqrt(2*pi) .* config.g_sigma)); % convert the distances to gaussian weighted values
            sbin = nansum(ds_norm); % the number of spikes which fall into the bin
            spikemap(yvec==ycen(bb),xvec==xcen(bb)) = sbin;        
        end
        ratemap = spikemap ./ dwellmap;
        config.dwellmap = dwellmap;
        config.spikemap = spikemap;     
        config.ratemap = ratemap; 
        
    case {3,'adaptive'}
        % this method was adapted from: Yartsev and Ulanovsky (2013) Representation of Three-Dimensional Space in the Hippocampus of Flying Bats
        % for every bin, we extend it outwards as a sphere until it contains at least config.min_dwell amount of time (or position data)
        % we then count how many spikes fall inside it and calculate firing rate as spikes/time. This method can be slow if there are many bins.
        % generally the ratemap will also have to be smoothed with a 3D gaussian. Although a dwellmap is technically produced it is not very
        % helpful as every bin will contain ~config.min_dwell amount of time.
     
        % prepare empty arrays
        ratemap = NaN(numel(yvec),numel(xvec));   
        if isfield(config,'dwellmap') && ~isempty(config.dwellmap)
            dwellmap = config.dwellmap;
        else        
            dwellmap = zeros(numel(yvec),numel(xvec));
        end
        spikemap = zeros(numel(yvec),numel(xvec));
        
        bpoints = ceil(config.min_dwell .* config.srate); % the number of position data points we need in a bin to fulfill binmin
        for bb = 1:length(xcen) % for every bin         
            if ~isfield(config,'dwellmap')
                d = sqrt((pox-xcen(bb)).^2 + (poy-ycen(bb)).^2); % get the distance to every position point
                ds = sort(d(:),'ascend');                
                if nanmin(ds) > config.min_dist % if this bin is further from the position data than our minimum cutoff, skip it
                    continue
                end
                mind = ds(bpoints+1);
                tbin = sum(d < mind)*(1/config.srate); % time in bin
                dwellmap(yvec==ycen(bb),xvec==xcen(bb)) = tbin;    
            else
                tbin = dwellmap(yvec==ycen(bb),xvec==xcen(bb));
            end
            
            if ~logical(tbin)
                continue
            end
            
            if ~numel(spx)
                ratemap(yvec==ycen(bb),xvec==xcen(bb)) = 0;
                continue
            end

            % calculate firing rate
            d2 = sqrt((spx-xcen(bb)).^2 + (spy-ycen(bb)).^2); % get the distance to every spike
            sbin = sum(d2 < mind); % the number of spikes which fall into the bin
            fbin = sbin ./ tbin; % the bin firing rate
            ratemap(yvec==ycen(bb),xvec==xcen(bb)) = fbin;     
            spikemap(yvec==ycen(bb),xvec==xcen(bb)) = sbin;                    
        end   
        config.ratemap = ratemap; 
        config.dwellmap = dwellmap; 
        config.spikemap = spikemap; 
        
    case {4, 'KDE'}
        % this method is entirely new. it uses a multivariate kernel density estimate approach to estimate probability density functions
        % for the position data and spike data separately. These are then converted to spike or dwelltime maps by multiplying by the total
        % of their respective variables. The firing rate map is then calculated as the spikemap divided by the dwelltime map. This ratemap
        % shouldn't need any smoothing, all smoothing should be achieved using the bandwidth of the KDE which can be adjusted. Larger bandwidth
        % is analogous to larger smoothing.
        f = mvksdensity([pox,poy],[xcen,ycen],'bandwidth',config.map_sigma,'Function','pdf','Kernel','normal');
        %f = akde([pox,poy],[xcen,ycen],50); % another ksdensity function that can estimate optimal bandwidth
        f = f .* config.bin_size; % convert to probability
        f = f .* (numel(pox).*(1/config.srate)); % convert to time        
        dwellmap = NaN(numel(yvec),numel(xvec));
        [~,idx] = ismember(xcen,xvec);
        [~,idy] = ismember(ycen,yvec);
        ida = sub2ind(size(dwellmap),idy,idx);
        dwellmap(ida) = f;
        dwellmap(dwellmap < config.min_dwell) = NaN;
        
        f = mvksdensity([spx,spy],[xcen,ycen],'bandwidth',config.map_sigma,'Function','pdf','Kernel','normal');        
        %f = akde([spx,spy],[xcen,ycen],50); % another ksdensity function that can estimate optimal bandwidth       
        f = f .* config.bin_size; % convert to probability
        f = f .* numel(spx); % convert to spikes           
        spikemap = NaN(numel(yvec),numel(xvec));
        ida = sub2ind(size(spikemap),idy,idx);
        spikemap(ida) = f;
        spikemap(isnan(dwellmap)) = NaN;
        
        ratemap = spikemap ./ dwellmap;
        ratemap(isnan(dwellmap)) = NaN;        
        config.dwellmap = dwellmap;
        config.spikemap = spikemap;     
        config.ratemap = ratemap; 
end

% adjust position and spike data to match maps
% config.poxnew = (pox - mean([min(pox) max(pox)])) ./ bin_size + (size(ratemap,2)/2);
% config.poynew = (poy - mean([min(poy) max(poy)])) ./ bin_size + (size(ratemap,1)/2);
% config.spxnew = (spx - mean([min(pox) max(pox)])) ./ bin_size + (size(ratemap,2)/2);
% config.spynew = (spy - mean([min(poy) max(poy)])) ./ bin_size + (size(ratemap,1)/2);








