function [pox,poy,poz,pot] = getLATTICEpath(snames,pos,config3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function just loads all of the 3D trajectory info produced by lattice for a session
%   given a set of filenames, an mtint (to check position data lengths) and some settings (optional)
%   [pox,poy,poz] = getLATTICEpath(snames,mtint,config3)
%
%%%%%%%% Inputs
%   snames = cell array of file names, no extensions
%   mtint = an mtint structure formed by klustest
%   config 3 = (optional) settings structure given by klustest3, otherwise defaults are used
%
%%%%%%%% Outputs
%   pox,poy,poz,pot = position x, y, z, t respectively
%
%%%%%%%% Comments
%   10/11/17 created so klustest can take advantage of reconstructed trajectories
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('config3','var') || isempty(config3) || all(isnan(config3(:)))
    % manual settings if config3 is missing
    samp_rate_hz = 50; % synchronised dacqTrack data should be 50Hz
    pos_tb = 1 / samp_rate_hz;    
    pixel_ratio = 1000;
    leds = 1; % typically we only use 1 LED
else
    % settings from config3 (given by klustest3)
    samp_rate_hz = config3.srate; % synchronised dacqTrack data should be 50Hz
    pos_tb = 1 / samp_rate_hz;    
    pixel_ratio = config3.pratio;
    leds = config3.leds; % typically we only use 1 LED
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the data
pox = []; 
poy = []; 
poz = []; % initialise variables
nsess = length(snames);
for pp = 1:nsess % for each session
    snow = snames{pp}; % get its filename
    fcheck = [pwd '\lattice\' snow '_merge_tracking.mat']; % load its 3d reconstructed path
    if exist(fcheck,'file') % check it exists
        load(fcheck);
    else
        error('Predetermined 3D trajectory not found: %s\n ...exiting',fcheck)
    end

    posn = [poxw(:) poyw(:) pozw(:)]; % extract the smoothed, weighted mean trajectory
    scount_now = pos.trial_samples(pp); % only take the first n samples, where n = the number of samples in the dacqUSB data (often there is one extra data point in the synchronised data)
    dcount = abs(length(posn(:,1))-scount_now);
    if length(posn(:,1)) < scount_now
        disp(sprintf('\t...3D session shorter by %d points (%.2fs), NaN padding data',dcount,(dcount * (1/samp_rate_hz))));
        posn2 = NaN(scount_now,3);
        posn2(1:size(posn,1),1:size(posn,2)) = posn;
        posn = posn2;
    elseif length(posn(:,1)) > scount_now
        disp(sprintf('\t...3D session longer by %d points (%.2fs), trimming data',dcount,(dcount * (1/samp_rate_hz))));   
        posn = posn(1:scount_now,:);        
    end

    pox = [pox; posn(:,1)]; % extract the x data and concatenate it
    poy = [poy; posn(:,2)]; % extract the y data and concatenate it
    poz = [poz; posn(:,3)]; % extract the z data and concatenate it
end
pox = double(pox(:));
poy = double(poy(:));
poz = double(poz(:));

