function [ratemap,dwellmap,spikemap,config,mapset] = mapDATA(pos,spk,config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION  function for mapping data such as position and spikes
% This function bins spike and position data into a firing rate map and can use various approaches
% This version has been reworked so that x,y data can be easily corrected to match bins in the map. 
% Data should be in cm before submitting to this function. Position data can be converted from pixels 
% to cm like this:
% position = position ./ pixels_per_metre .* 100;
%
% USAGE:
%         [ratemap,config] = mapDATA5(pos,spk,...)
%
% INPUT:
%         pos - position data to make map with [x,y]
%         spk - spike data to make map with [x,y]
%         config - input configuration structure containing the optional fields:
%             dwellmap - if provided this will be used in place of generating a new one for a speed increase
%             rmethod - ratemap method, choices include 'nearest', 'gaussian', 'adaptive' or 'KDE', default is 'nearest'
%             smethod - smoothing method, 1 = smooth(spk/time), 2 = smooth(spk)/smooth(time), 3 = no smoothing at all, used by nearest method only, default is 1
%             map_padd - how many bins to padd maps with, default is 2
%             map_lims - the limits of the map, [min_x max_x min_y max_y], default is the limits of the position data
%             bin_size - the bin size in cm to use, default is 2
%             map_sigma - the sigma of the gaussian kernel used to smooth the maps, default is 2
%             min_dwell - the minimum time that must be in a bin for it to be considered visited, default is 0.1
%             g_sigma - sigma for the gaussian weighting of position/spike data, used by gaussian method only, default is 10
%             min_dist - minimum distance bin has to be from position data for it to be considered visited, default is 1
%             max_dist - the maximum distance data can be from a bin and still influence its value, used by gaussian method only, default is 3
%             srate - the sampling rate of the data in Hz, used to calculate time, default is 50
%
% OUTPUT:
%    ratemap - generated firing rate map (Hz)
%    dwellmap - generated firing rate map (seconds)
%    spikemap - generated firing rate map (spikes)
%    output - structure containing additional data and the settings used to generate the map
%
% EXAMPLES:
%     % generate a ratemap using default settings
%     pos = rand(500000,2)*200; % positions in cm (a 2 m environment)
%     spk = [normrnd(100,10,5000,1) normrnd(100,20,5000,1)]; % spikes in cm
%     [r,d,s,o] = mapDATA(pos,spk);
%     figure
%     subplot(2,2,1)
%     plot(pos(:,1),pos(:,2),'k'); hold on;
%     plot(spk(:,1),spk(:,2),'r.')
%     daspect([1 1 1])
%     subplot(2,2,2)
%     imagesc(r);
%     daspect([1 1 1])
% 
%     % use the dwellmap from the above
%     config = struct;
%     config.dwellmap = d;
%     [r,d,s,o] = mapDATA(pos,spk,config);
%     subplot(2,2,3)
%     imagesc(r);
%     daspect([1 1 1])
% 
%     % make a new map with different settings
%     config = struct;
%     config.bin_size = 5;
%     config.map_sigma = 3;
%     [r,d,s,o] = mapDATA(pos,spk,config);
%     subplot(2,2,4)
%     imagesc(r);
%     daspect([1 1 1])
%
% See also: KLUSTEST HIST3 IMGAUSSFILT

% HISTORY:
% version 1.0.0, Release 08/08/16 created as a simple mapping function
% version 2.0.0, Release 20/04/17 created from mapDATA because I want a function I can use with overdispersion
% version 2.0.1, Release 20/04/17 function now outputs converted position data (to match map)
% version 2.0.2, Release 20/04/17 added limit capability so spike and time maps can be made with same limits
% version 3.0.0, Release 27/07/17 changed to v3, this version expects data to be in cm, not pixels (meant to work with variable pixel ratio klustest)
% version 4.0.0, Release 01/08/17 changed to v4, overhauled, added switches for different approaches
% version 4.1.0, Release 02/08/17 added adaptive method
% version 4.2.0, Release 03/08/17 added KDE approach, comments
% version 4.2.1, Release 31/10/17 added better functionality for converting units of data
% version 5.0.0, Release 08/08/18 overhauled for klustest update, simplified inputs, removed redundant code
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
field_names = {'dwellmap','rmethod','smethod','map_padd','map_lims','bin_size','map_sigma','min_dwell','g_sigma','min_dist','max_dist','srate'};
default_vals = {[],'nearest',1,2,[],2,1.5,0.1,10,1,3,50};
for ff = 1:length(field_names)
    if ~isfield(config,field_names{ff})
        config.(field_names{ff}) = default_vals{ff};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
% extract data
pox = single(pos(:,1));
poy = single(pos(:,2)); 
spx = single(spk(:,1));
spy = single(spk(:,2)); 

% prepare limits of map
if isempty(config.map_lims) || all(isnan(config.map_lims))
    lx = [nanmin(pox) nanmax(pox)];
    ly = [nanmin(poy) nanmax(poy)];
else
    lx = config.map_lims(1:2);
    ly = config.map_lims(3:4);    
end
mapset = struct;
mapset.map_lims = [lx ly];

%% Prepare bins
% centre data on the origin
pox = pox - mean(lx); 
poy = poy - mean(ly);   
spx = spx - mean(lx); 
spy = spy - mean(ly);   
lx = lx - mean(lx); 
ly = ly - mean(ly);  

% vectors for bin edges (used by nearest, gaussian, adaptive methods)
xvec = 0 : config.bin_size : (max(abs(lx))+config.bin_size+config.bin_size*config.map_padd); 
yvec = 0 : config.bin_size : (max(abs(ly))+config.bin_size+config.bin_size*config.map_padd);    

% sort them to increase from -max to max
xvec = unique(sort([-xvec xvec],'ascend')); 
yvec = unique(sort([-yvec yvec],'ascend'));
mapset.xloc = xvec;
mapset.yloc = yvec;

mapset.poxnew = single((pox ./ config.bin_size) + ceil(length(xvec)/2));
mapset.poynew = single((poy ./ config.bin_size) + ceil(length(yvec)/2));
mapset.spxnew = single((spx ./ config.bin_size) + ceil(length(xvec)/2));
mapset.spynew = single((spy ./ config.bin_size) + ceil(length(yvec)/2));

if isfield(config,'data_to_process')
    for i = 1:length(config.data_to_process)
        datnow = config.data_to_process{i};
        xnow = datnow(:,1);
        ynow = datnow(:,2);
        xnow = xnow - mean(mapset.map_lims(1:2)); 
        ynow = ynow - mean(mapset.map_lims(3:4));  
        xnow = single((xnow ./ config.bin_size) + ceil(length(xvec)/2));
        ynow = single((ynow ./ config.bin_size) + ceil(length(yvec)/2));
        mapset.dataout{i} = [xnow(:) ynow(:)];
    end
end

% determine if we will need a new dwellmap or not
gen_dwell = 1; % we are going to have to make a new dwellmap
if isfield(config,'dwellmap')
    dwellmap = config.dwellmap;    
    gen_dwell = 0; % we don't have to make a new dwellmap
    
    if isempty(dwellmap) || ~all(isnan(dwellmap(:)))
        gen_dwell = 1; % we are going to have to make a new dwellmap
    end
end
rmap = NaN(size(dwellmap));
smap = zeros(size(dwellmap));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE MAP
switch config.rmethod
    case {1,'nearest'} % 3D histogram
        % this method is the most widely used and can be found everywhere in the literature. It is also the most basic and the fastest.
        % data are literally binned where they are found. This is done for position data and the bin counts are then multiplied by the
        % sampling rate to get time. The same is then done for spikes. The spikemap and dwellmap are then gaussian smoothed and then
        % the spikemap is divided by the dwellmap. Empty bins, or bins where the rat spent less than config.min_dwell time are empty (NaN)
        
        % generate spikemap
        spikemap = hist3([spy,spx],{yvec,xvec});

        % generate dwellmap
        if gen_dwell
            dwellmap = hist3([poy,pox],{yvec,xvec});             
        end
        
        % generate and smooth ratemap
        if config.smethod==1 && config.map_sigma>0
            spikemap = imgaussfilt(spikemap,config.map_sigma);
            dwellmap = imgaussfilt(dwellmap,config.map_sigma) .* (1/config.srate);
            ratemap = spikemap ./ dwellmap;          
        elseif config.smethod==2 && config.map_sigma>0
            dwellmap = dwellmap .* (1/config.srate);
            ratemap = spikemap ./ dwellmap;
            ratemap = imgaussfilt(ratemap,config.map_sigma); 
        elseif config.smethod==3 || config.map_sigma==0
            dwellmap = dwellmap .* (1/config.srate);
            ratemap = spikemap ./ dwellmap;            
        end

        % make sure there are no weird values in the ratemap
        ratemap(dwellmap < config.min_dwell | dwellmap == 0 | isnan(dwellmap) | isinf(dwellmap)) = NaN;
        
    case {2,'gaussian'}
        % this is an approach taken from: Leutgeb et al. (2005) Independent Codes for Spatial and Episodic Memory in Hippocampal Neuronal Ensembles
        % for every bin we calculate the distance to every spike and every position data point. These distances are weighted with a gaussian so that
        % data close to the bin centre = 1 and data far away = 0. The sum of the spike distances are then divided by the sum of the position
        % distances to obtain firing rate. The standard deviation of the gaussian can be changed, the result being that bins become 'larger' or rather
        % bins will overlap more as more data is contained in each. This method can be slow when there are many bins.
        
        % prepare bin centres
        [X1,Y1] = meshgrid(xvec,yvec); 
        xcen = X1(:);
        ycen = Y1(:);
        mapset.xv = xcen;
        mapset.yv = ycen;

        ratemap = NaN(numel(yvec),numel(xvec));           
        if gen_dwell
            dwellmap = zeros(numel(yvec),numel(xvec));
        end
        spikemap = zeros(numel(yvec),numel(xvec));

        for bb = 1:length(xcen) % for every bin   
            if gen_dwell
                % get the distance to every position point
                rindx = pox>xcen(bb)-config.max_dist & pox<xcen(bb)+config.max_dist & poy>ycen(bb)-config.max_dist & poy<ycen(bb)+config.max_dist;            
                dp = sqrt(sum(([pox(rindx) poy(rindx)]-[xcen(bb) ycen(bb)]).^2,2));
                if nanmin(dp) > config.min_dist % if this bin is too far from the position data
                    continue
                end 
                dp_norm = (exp(-0.5 * ((dp - 0)./config.g_sigma).^2) ./ (sqrt(2*pi) .* config.g_sigma)) ./ (1 ./ (sqrt(2*pi) .* config.g_sigma)); % convert the distances to gaussian weighted values
                tbin = nansum(dp_norm)*(1/config.srate); % time in voxel
                dwellmap(yvec==ycen(bb),xvec==xcen(bb)) = tbin; % accumulate data
            else
                tbin = dwellmap(yvec==ycen(bb),xvec==xcen(bb));
            end

            % get the distance to every spike
            if tbin < config.min_dwell % if the rat spent too little time in this voxel
                spikemap(yvec==ycen(bb),xvec==xcen(bb)) = 0;
                continue
            end
            rindx = spx>xcen(bb)-config.max_dist & spx<xcen(bb)+config.max_dist & spy>ycen(bb)-config.max_dist & spy<ycen(bb)+config.max_dist;            
            ds = sqrt(sum(([spx(rindx) spy(rindx)]-[xcen(bb) ycen(bb)]).^2,2));
            ds_norm = (exp(-0.5 * ((ds - 0)./config.g_sigma).^2) ./ (sqrt(2*pi) .* config.g_sigma)) ./ (1 ./ (sqrt(2*pi) .* config.g_sigma)); % convert the distances to gaussian weighted values
            sbin = nansum(ds_norm); % the number of spikes which fall into the bin
            spikemap(yvec==ycen(bb),xvec==xcen(bb)) = sbin;        
        end
        ratemap = spikemap ./ dwellmap;

        % make sure there are no weird values in the ratemap
        ratemap(dwellmap < config.min_dwell | dwellmap == 0 | isnan(dwellmap) | isinf(dwellmap)) = NaN;        
        
    case {3,'adaptive'}
        % this method was adapted from: Yartsev and Ulanovsky (2013) Representation of Three-Dimensional Space in the Hippocampus of Flying Bats
        % for every bin, we extend it outwards as a sphere until it contains at least config.min_dwell amount of time (or position data)
        % we then count how many spikes fall inside it and calculate firing rate as spikes/time. This method can be slow if there are many bins.
        % generally the ratemap will also have to be smoothed with a 3D gaussian. Although a dwellmap is technically produced it is not very
        % helpful as every bin will contain ~config.min_dwell amount of time.
     
        % prepare bin centres
        [X1,Y1] = meshgrid(xvec,yvec); 
        xcen = X1(:);
        ycen = Y1(:);
        mapset.xv = xcen;
        mapset.yv = ycen;
        
        % prepare empty arrays
        ratemap = NaN(numel(yvec),numel(xvec));   
        if gen_dwell      
            dwellmap = zeros(numel(yvec),numel(xvec));
        end
        spikemap = zeros(numel(yvec),numel(xvec));
        
        bpoints = ceil(min_dwell .* config.srate); % the number of position data points we need in a bin to fulfill binmin
        for bb = 1:length(xcen) % for every bin         
            if gen_dwell
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
        
        % make sure there are no weird values in the ratemap
        ratemap(dwellmap < config.min_dwell | dwellmap == 0 | isnan(dwellmap) | isinf(dwellmap)) = NaN;         
        
    case {4, 'KDE'}
        % this method is entirely new. it uses a multivariate kernel density estimate approach to estimate probability density functions
        % for the position data and spike data separately. These are then converted to spike or dwelltime maps by multiplying by the total
        % of their respective variables. The firing rate map is then calculated as the spikemap divided by the dwelltime map. This ratemap
        % shouldn't need any smoothing, all smoothing should be achieved using the bandwidth of the KDE which can be adjusted. Larger bandwidth
        % is analogous to larger smoothing.
        
        % prepare bin centres
        [X1,Y1] = meshgrid(xvec,yvec); 
        xcen = X1(:);
        ycen = Y1(:);
        mapset.xv = xcen;
        mapset.yv = ycen;
        
        % generate the dwell time map
        if gen_dwell     
            dwellmap = NaN(numel(yvec),numel(xvec));
            
            f = mvksdensity([pox,poy],[xcen,ycen],'bandwidth',config.map_sigma,'Function','pdf','Kernel','normal');
            f = f .* config.bin_size; % convert to probability
            f = f .* (numel(pox).*(1/config.srate)); % convert to time        
            [~,idx] = ismember(xcen,xvec);
            [~,idy] = ismember(ycen,yvec);
            ida = sub2ind(size(dwellmap),idy,idx);
            dwellmap(ida) = f;
            dwellmap(dwellmap < config.min_dwell) = NaN;
        end
        
        % generate the spike map
        spikemap = NaN(numel(yvec),numel(xvec));        
        f = mvksdensity([spx,spy],[xcen,ycen],'bandwidth',config.map_sigma,'Function','pdf','Kernel','normal');        
        f = f .* config.bin_size; % convert to probability
        f = f .* numel(spx); % convert to spikes           
        ida = sub2ind(size(spikemap),idy,idx);
        spikemap(ida) = f;
        spikemap(isnan(dwellmap)) = NaN;
        
        ratemap = spikemap ./ dwellmap;
        ratemap(isnan(dwellmap)) = NaN;        

        % make sure there are no weird values in the ratemap
        ratemap(dwellmap < config.min_dwell | dwellmap == 0 | isnan(dwellmap) | isinf(dwellmap)) = NaN;        
        
end





