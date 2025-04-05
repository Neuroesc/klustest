function sdata = getSDATA(mtint,tet,clu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an mtint file, performs a number of analyses and outputs an sdata structure
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial setup
clear sdata
sdata = struct;

% get initial data from mtint
clu_s = ['c_' num2str(clu)]; % string for use in structure array

spik_count = numel(mtint.tetrode(tet).pos_sample); % retrieve the number of spikes recorded on this tetrode
clus_count = numel(unique(mtint.tetrode(tet).cut));
sdata.(tet_s).clus_count = clus_count; % add data to structure
sdata.(tet_s).spik_count = spik_count; % add data to structure
   
% get some vectors that we can use to sort data
clu_identity = mtint.tetrode(tet).cut; % clu_assign is a vector of numbers, one for each spike, each number corresponds to a cluster
clu_indx = (clu_identity == clu); % logical vector, one number per spike, 1 if in cluster, 0 if not
pos_identity = mtint.tetrode(tet).pos_sample; % pos_assign is a vector of numbers, one for each spike, each number corresponds to a position data point
n_spikes = length(find(clu_identity == clu));

% calculate cluster firing rate
frate = n_spikes ./ duration;
sdata.(tet_s).(clu_s).frate = frate; % add data to structure

% retrieve position data
pox = sdata.pox; % add data to structure
poy = sdata.poy; % add data to structure
pot = sdata.pot; % add data to structure
head_direct = sdata.hd; % add data to structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n_spikes > 0 % if there is at least one spike
    %% Waveforms
    waves{1} = mtint.tetrode(tet).ch1;
    waves{2} = mtint.tetrode(tet).ch2;
    waves{3} = mtint.tetrode(tet).ch3;
    waves{4} = mtint.tetrode(tet).ch4;

    mean_wav = {}; % will hold mean waveform data
    std_wav = {}; % will hold standard deviation of waveform data
    max_wav = repmat(NaN,1,4); % will hold maximum amps of waveforms
    min_wav = repmat(NaN,1,4); % will hold minimum amps of waveforms
    width_wav = repmat(NaN,1,4); % will hold waveform widths
    for w = 1:4 % for every recording channel
        wav = waves{w};
        waves{w} = double(wav(clu_identity == clu,:));

        ch = squeeze(mean(waves{w},1));
        chs = squeeze(std(waves{w},1));
        [maxval,maxt] = max(ch);
        [postminval postmint] = min(ch(maxt:end));
        postmint = postmint + maxt - 1;
        width = postmint - maxt;
        width = width * (1000/50); 
        width_wav(w) = width;
        max_wav(w) = maxval;
        min_wav(w) = postminval;
        mean_wav{w} = ch;
        std_wav{w} = chs;
    end % for w = 1:4 % for every recording channel

    % accumulate data
    sdata.(tet_s).(clu_s).wave_widths = width_wav; % add data to structure
    sdata.(tet_s).(clu_s).wave_maxs = max_wav; % add data to structure
    sdata.(tet_s).(clu_s).wave.mins = min_wav; % add data to structure
    sdata.(tet_s).(clu_s).wave_stds = std_wav; % add data to structure


    %% Spike data
    spx = pox(pos_identity(clu_indx)); % the x coordinate of every spike in this cluster
    spy = poy(pos_identity(clu_indx)); % the y coordinate of every spike in this cluster
    spt = spiketime(clu_indx); % the time point for every spike in this cluster

    % sort out head direction data
    hd = double(head_direct(pos_identity(clu_identity == clu))); % get an index of which spikes belong to this cluster, then get an indext of whcih positions these spikes correspond to, then use this to get the correct head directions       

    % accumulate data
    sdata.(tet_s).(clu_s).spike_spx = spx; % add data to structure
    sdata.(tet_s).(clu_s).spike_spy = spy; % add data to structure
    sdata.(tet_s).(clu_s).spike_spt = spt; % add data to structure
    sdata.(tet_s).(clu_s).hd = hd; % add data to structure


    %% Firing rate map and analyses
    map_limits = [min(pox)-map_padd max(pox)+map_padd min(poy)-map_padd max(poy)+map_padd];

    % create dwell time map
    [dwellmap,~] = mapDATA(pox,poy,map_limits,bin_size,pixel_ratio);
    dwellmap = dwellmap .* pos_tb;
    dwellmap(dwellmap < min_dwell) = 0;
    dwellmap = imgaussfilt(dwellmap,map_sigma);
    sdata.(tet_s).(clu_s).spatial_dwellmap = dwellmap; % add data to structure

    % create spike map
    [spikemap,~] = mapDATA(spx,spy,map_limits,bin_size,pixel_ratio);
    spikemap = imgaussfilt(spikemap,map_sigma);

    % create ratemap
    ratemap = spikemap ./ dwellmap;
    ratemap(dwellmap == 0) = NaN;

    % ratemap analysis
    skaggs = skaggs_info2(ratemap,dwellmap);
    spars = sparsity(ratemap,dwellmap);
    cohe = spatial_coherence(ratemap,dwellmap);

    % place field analysis
    [fieldd] = getPFIELDS(ratemap);

    % accumulate data
    sdata.(tet_s).(clu_s).spatial_ratemap = ratemap; % add data to structure
    sdata.(tet_s).(clu_s).spatial_information = skaggs; % add data to structure
    sdata.(tet_s).(clu_s).spatial_sparsity = spars; % add data to structure
    sdata.(tet_s).(clu_s).spatial_coherance = cohe; % add data to structure
    sdata.(tet_s).(clu_s).field_count = length(fieldd.fields(:,1)); % add data to structure
    sdata.(tet_s).(clu_s).field_data = fieldd; % add data to structure
              
    %% Grid autocorrelation and analyses
    % create autocorrelation
    automap = GridAutoCorr(ratemap);

    % autocorrelation analysis
    [grid_score,grid_spacing,field_size,grid_orientation,grid_ellipticity] = GridAnalysis(automap,bin_size);

    % accumulate data
    sdata.(tet_s).(clu_s).grid_autocorrelation = automap; % add data to structure
    sdata.(tet_s).(clu_s).grid_score = grid_score; % add data to structure
    sdata.(tet_s).(clu_s).grid_spacing = grid_spacing; % add data to structure
    sdata.(tet_s).(clu_s).grid_field_size = field_size; % add data to structure
    sdata.(tet_s).(clu_s).grid_orientation = grid_orientation; % add data to structure
    sdata.(tet_s).(clu_s).grid_ellipticity = grid_ellipticity; % add data to structure

    %% Head direction analysis
    ai = linspace(0,2*pi,hd_bins)'; % angles for binning
    if strcmp(hd_type,'density')
        hd_s = deg2rad(head_direct); % session head direction in radians
        hd_c = deg2rad(hd); % cell head direction in radians
        [hd1] = circ_ksdensity(hd_s,ai,[],0.02); % the session head direction       
        [hd2] = circ_ksdensity(hd_c,ai,[],0.02); % the cell's head direction
        hd1 = hd1 .* pos_tb; % convert session HD to time
        hd3 = hd2 ./ hd1; % calculate HD firing rate

    elseif strcmp(hd_type,'histogram')
        hd1 = hist(deg2rad(hd),hd_bins); % the session head direction   
        hd2 = hist(deg2rad(head_direct),hd_bins); % the cell's head direction
        % nAll(nAll < sample_rate) = nan; % min dwell cutoff
        hd1 = hd1 .* pos_tb; % convert session HD to time
        hd3 = hd2 ./ hd1; % calculate HD firing rate
        fh = fspecial('average',[1 5]);
        hd3 = imfilter(hd3,fh,'circular','same');

    end % if strcmp(hd_type,'density')
    hd1 = hd1 ./ max(hd1); % normalise session hd
    hd3 = hd3 ./ max(hd3); % normalise cell hd

    % head direction analyses
    rayleigh = circ_r(ai,hd3); % rayleigh vector length
    mx2 = rad2deg(ai(hd3 == max(hd3))); % preferred angle (location of max frate)

    % accumulate data
    sdata.(tet_s).(clu_s).hd_frate = hd_c; % add data to structure
    sdata.(tet_s).(clu_s).hd_rayleigh = rayleigh; % add data to structure
    sdata.(tet_s).(clu_s).hd_maximum = mx2; % add data to structure

    %% Spikes vs time 
    tvals = (0:time_bins:duration); % vector of time points at which we should calculate spike probability
    [c,b] = histc(spt,tvals);
    [probs,xi] = ksdensity(spt,tvals);

    % accumulate data
    sdata.(tet_s).(clu_s).spikes_time_histogram = [c,b]; % add data to structure
    sdata.(tet_s).(clu_s).spikes_time_ksdensity = [probs(:),xi(:)]; % add data to structure

    %% Spike phase analysis
    % phase preference plot
    lfpf = mtint.lfp(1).theta;
    sdata.theta = lfpf; % add data to structure    
    t = mtint.lfp.t;

    % estimate theta phase, based on lfp data
    [phase_out,~,~] = Phase([t lfpf],spt);
    spp = phase_out(:,2);

    % bin the theta phase data
    ai = -pi:0.1:pi;
    yi = histc(spp,ai);
    [probs,xi] = ksdensity([spp; (spp+2*pi)],ai);

    % example cosine wave
    swav = cos(ai);
    swav = ((swav ./ max(swav))+abs(min((swav ./ max(swav))))).*max(yi);

    mu = circ_mean(spp);

    % accumulate data
    sdata.(tet_s).(clu_s).phase_mean = mu; % add data to structure
    sdata.(tet_s).(clu_s).spike_phase = spp; % add data to structure
    sdata.(tet_s).(clu_s).spike_phase_ideal = swav; % add data to structure
    sdata.(tet_s).(clu_s).spike_phase_ksdensity = [probs,xi]; % add data to structure

    %% Phase map and analyses               
    if strcmp(pmap_type,'tri')
        warning('off','MATLAB:delaunay:DupPtsDelaunayWarnId')
        tri = delaunay(spx,spy);
        phmap = tri;
    elseif strcmp(pmap_type,'histogram')
        phasemap = mapDATA3(spx,spy,spp,map_padd,pmap_bs,pixel_ratio);
        nindx = find(isnan(phasemap));
        phasemap(nindx) = 0;
        phasemap = imgaussfilt(phasemap,map_sigma);
        phasemap(nindx) = NaN;
        phasemap = phasemap';
        phmap = phasemap;
    end % if strcmp(pmap,'tri')

    % accumulate data
    sdata.(tet_s).(clu_s).phase_map = phmap; % add data to structure

    %% Spike autocorrelation - theta analaysis           
    [corrdata1,tms] = spk_acorr2(spt,1,500);

    % fit a decomposing sine wave to the spike autocorrelogram
    [thetaR,thetaP,thetaIndx,thetaPowr,thetaLin] = getTHETAfit([tms' corrdata1]);

    % accumulate data
    sdata.(tet_s).(clu_s).theta_r = thetaR; % add data to structure
    sdata.(tet_s).(clu_s).theta_p_value = thetaP; % add data to structure
    sdata.(tet_s).(clu_s).theta_index = thetaIndx; % add data to structure
    sdata.(tet_s).(clu_s).theta_power = thetaPowr; % add data to structure
          
    %% Spike autocorrelation - refractory period analysis
    [corrdata2,tms] = spk_acorr2(spt,1,50);

    % Calculate refractory contamination
    % Navratilova and McNaughton (2016) Grids from bands, or bands from grids? An examination of the effects of single unit contamination on grid cell firing fields
    % Fee, Mitra, Kleinfeld (1996) Automatic sorting of multiple unit neuronal signals in the presence of anisotropic and non-Gaussian variability
    half_spike = corrdata2(tms >= 0); % take only the positive side of spike autocorrelogram
    half_time = tms(tms >= 0); % take only the positive side of spike autocorrelogram
    tau_tot = tau_r - tau_c;
    Ns = n_spikes; % number of cluster spikes
    lambda = Ns/duration;  % mean firing rate for cluster 
    RPV = sum((half_spike ~= 0) & (half_time' <= tau_r)); % get the number of refractory period violations

    % get Poisson confidence interval on number of expected RPVs
    conf_int = 95; % percent confidence interval
    [~,interval] = poissfit(RPV,(100-conf_int)/100); 

    % convert contamination from number of RPVs to a percentage of spikes
    RPVs = [RPV interval(1) interval(2)];
    cont_bounds = NaN(1,3);
    for i = 1:length(RPVs)
        RPVnow = RPVs(i);
        RPVT = 2 * tau_tot * Ns; % total amount of time in which an RPV could occur = the usable refractory period (tau_tot) around each spike (*Ns) and on each side of the spike (*2)
        RPV_lambda = RPV / RPVT; % rate of RPV occurence
        p =  RPV_lambda / lambda; % estimate of % contamination of cluster

        % force p to be a real number in [0 1]
        if isnan(p)
            p = 0; 
        end % if isnan(p)
        if p > 1
            p = 1; 
        end   % if p > 1
        cont_bounds(i) = p;
    end % for i = 1:length(RPVs)
    censored_estimate = (tau_c/1000) * spik_count / duration; % estimate the fraction of spikes not detected because of the system lockout

    % accumulate data
    sdata.(tet_s).(clu_s).refractory_violations = RPV; % add data to structure
    sdata.(tet_s).(clu_s).refractory_contamination = cont_bounds(1); % add data to structure
    sdata.(tet_s).(clu_s).refractory_contamination_95_lower_upper_bounds = cont_bounds(2:3); % add data to structure
    sdata.(tet_s).(clu_s).refractory_censoring = cont_bounds(2:3); % add data to structure

    %% Mahalanobis distance, cluster quality analyses
    [isod,lratio,nd,cd,fdata,nfets] = clusterQUALITY(cname,tet,clu);

    % accumulate data
    sdata.(tet_s).(clu_s).cluster_isolation = isod; % add data to structure
    sdata.(tet_s).(clu_s).cluster_lratio = lratio; % add data to structure                    

    %% Cell type/identity
    if ignore_pos
        frate = sdata.(tet_s).(clu_s).frate;
        max_wav = sdata.(tet_s).(clu_s).wave_maxs;
        max_wide = sdata.(tet_s).(clu_s).wave_widths;
        [v,windx] = max(max_wav);
        width = max_wide(windx);

        c_type = [];
        if width > 250 % putative pyramidal cells
            if frate > 5
                c_type = 'high pyramidal';
            elseif frate < 0.1
                c_type = 'silent cell';
            else
                c_type = 'medium pyramidal';
            end % if frate > 5
        elseif width < 250 % putative interneurons
            c_type = 'Interneuron';
        end % if width > 250 % putative pyramidal cells 

    else
        frate = sdata.(tet_s).(clu_s).frate;
        max_wav = sdata.(tet_s).(clu_s).wave_maxs;
        max_wide = sdata.(tet_s).(clu_s).wave_widths;
        skaggs = sdata.(tet_s).(clu_s).spatial_information;
        gscore = sdata.(tet_s).(clu_s).grid_score;
        rvect = sdata.(tet_s).(clu_s).hd_rayleigh;

        [v,windx] = max(max_wav);
        width = max_wide(windx);

        c_type = [];
        if width > 250 % putative pyramidal cells
            if frate > 5
                c_type = 'high pyramidal';
            elseif frate < 0.1
                c_type = 'silent cell';
            else
                if skaggs > 0.5
                    c_type = 'place cell';
                    c_count(1) = c_count(1) + 1;
                else
                    c_type = 'medium pyramidal';
                end % if skaggs > 0.5
            end % if frate > 5
        elseif width < 250 % putative interneurons
            c_type = 'Interneuron';
        end % if width > 250 % putative pyramidal cells

        if grid_score > 0.2
            c_type = [c_type ' & grid cell'];
            c_count(2) = c_count(2) + 1;
        end % if grid_score > 0.2

        if rayleigh > 0.2
            c_type = [c_type ' & hd cell'];
            c_count(3) = c_count(3) + 1;
        end % if grid_score > 0.2     
    end % if ignore_pos

    % accumulate data
    sdata.(tet_s).(clu_s).cell_type = c_type; % add data to structure

end % if n_spikes > 0 
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    



