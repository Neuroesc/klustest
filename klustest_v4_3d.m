% function kilotest(varargin)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  short description
% long description
%
% USAGE:
%       [out] = template(in,in2)
%
% INPUT:
%       in - input 1
%       in2 - input 2
%
% OUTPUT:
%       p - output
%
% EXAMPLES:
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    p = inputParser;
    addParameter(p,'tetrodes',1:16,@(x) ~isempty(x) && ~all(isnan(x(:))) && isnumeric(x));  
    addParameter(p,'clusters',0,@(x) ~isempty(x) && ~all(isnan(x(:))) && isnumeric(x)); 
    addParameter(p,'cname','kwiktint',@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    dlist = strsplit(pwd,'\');
    rname = dlist{end-1};
    addParameter(p,'rname',rname,@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    
%     parse(p,varargin{:});
    parse(p);    
    config = p.Results;  

    %% overrides
    override.part_config = 0;
    override.use_groups = 1;
    
    %% Map settings
    mapset.ppm          = 300;
    mapset.method       = 'histogram';
    mapset.binsize      = 20; % firing rate map bin size
    mapset.ssigma       = 40; % firing rate map smoothing sigma
    mapset.padding      = 10; % in mm
    mapset.mindwell     = 0.01; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    mapset.mindist      = 40; % (mm, default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    mapset.smethod      = 1; % smoothing method, 1 = before division, 2 = after, 3 = no smoothing    
    mapset.frcut        = 0.2; % cutoff % for regions to be considered a place field
    mapset.arcut        = 36; % cutoff area for regions to be considered a place field
    
    % HD settings - settings to use when generating head direction maps
    mapset.hd_type      = 'histogram'; % (default 'histogram') enter 'density' for a circular kolmogorov smirnov density estimate plot, enter 'histogram' for the traditional histogram polar plot
    mapset.hd_bins      = 64; % (default 64) the number of bins to use when computing HD plot
    mapset.hd_sigma     = 0.04; % (default 2) the standard deviation of the gaussian kernel used in HD circular density estimate
    mapset.hd_boxcar    = 3; % (defualt 3) the number of bins over which to compute the HD histogram boxcar average    
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PARTITION SETTINGS
    % part_names gives the names you want the outputs to be saved as, these cannot start with a number: i.e. 'session1'
    part_names = {'arena1','hills','arena2'};
    % part_methods specifies the method of partition, corresponding to each of the names above: 1 = combine everything, 2 = take a recording session, 3 = use digital inputs
    % part_methods must be a cell array with the same length as part_names    
    part_methods = {2 2 2};
    % part_intervals gives vector(s) indicating which intervals to include in each partition: 
    % if method = 1 this value is ignored
    % if method = 2 this should specify which recording sessions to include
    % if method = 3 this should specify which digital input pairs to include (inf = use all) 
    % part_intervals must be a cell array with the same length as part_names    
    % i.e. if part_methods={2; 2; 2}, then if part_intervals={1; 2; 3}, rec 1 will go in part 1, rec 2 in part 2 and rec 3 in part 3
    % if part_methods={2; 2; 2; 2; 2}, then if part_intervals={1; [2 4 5]; 3}, rec 1 will go in part 1, rec 2,4 and 5 in part 2 and rec 3 in part 3 
    % if part_methods={2; 3}, then if part_intervals={2; [1 2 3 5 6 7 12]}, rec 2 will go in part 1, part 2 will be made up of the time segments between interval key pairs 1 2 3 5 6 7 12
    part_intervals = {1 2 3};
    % part_interval_keys specifies what keypresses were used to delineate trial starts and ends {start,end} these can be integers or characters
    % i.e. {1,2} or {'s','e'}, {'1','2'} is the same as {1,2}
    % part_interval_keys must be a cell array with the same length as part_names        
    part_interval_keys = {{},{},{}};
    % convert this information into a convenient table (which we will add to later)
    part_config = table(part_names,part_methods,part_intervals,part_interval_keys);

    % The part_config table contains basic information about how/where/if to separate the data into different
    % partitions or parts. For instance if we record an open field, then some sort of maze, then the open field
    % again, we may want to divide the data into those 3 parts for separate analysis, or we might want to just 
    % lump everything together if we care about some intrinsic property of the cells
    % We want to save this table so that in the future we can just run the function with the same settings    
    [~,~,~] = mkdir([pwd '\' config.cname]); % create a folder to hold outputs        
    part_config_fname = [pwd '\' config.cname '\part_config'];
    if override.part_config || ~exist(part_config_fname,'file')
        writetable(part_config,part_config_fname)
    elseif override.part_config && exist(part_config_fname,'file')
        part_config_fname2 = [pwd '\' config.cname '\part_config_' datestr(now,30) '.mat'];      
        [~,~,~] = movefile(part_config_fname,part_config_fname2);
        writetable(part_config,part_config_fname)
    elseif ~override.part_config && exist(part_config_fname,'file')
        part_config = readtable(part_config_fname);
    end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PREPARE DATA
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Tetrodes and sessions
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'); tic;
    stk = dbstack;
    tnow = datestr(now,'yyyy-mm-dd-HH-MM-SS');
    disp(sprintf('Running %s at %s...',stk.name,tnow))

%% >>>>>>>>>> find which tetrodes are available
    % When sessions are cluster cut together (which you should do with multiple sessions recording the same cells)
    % kilocut saves the filenames of the individual sessions in the kilo.mat file as well as the tetrodes analysed
    disp(sprintf('Assessing data...'))
    [~,snames,cutname] = get_tets_for_klustest('Neuralynx',config); 
    nsess = size(snames,1);
    data_dirs = strcat(repmat(pwd,nsess,1),repmat('\',nsess,1),snames); % where the original cheetah files are stored
    disp(sprintf('\t...read %s',cutname));
    disp(sprintf('\t...%d sessions',nsess)); 
    
    % find which tetrodes we want to (and can) analyse
    possible_tetrodes = config.tetrodes;
    for ff = 1:length(data_dirs) % for every data directory
        d = dir([data_dirs{ff} '\*.ntt']);
        fnames = {d.name};
        tets_found = cell2mat(cellfun(@str2num,replace(fnames,{'TT','.ntt'},''),'UniformOutput',false));
        possible_tetrodes = intersect(possible_tetrodes,tets_found);
    end
    tets = possible_tetrodes; 
    disp(sprintf('\t...will analyse tetrodes: %s',mat2str(tets)))                                 

    % Start to prepare the pdata (part or session data) table
    % To save space and make things easier for the user I have divided the data into two main arrays
    % the pdata structure contains data which applies to the whole session, like positions or dwell maps
    % the sdata table contains data for each cluster/part (one per row). The idea being that tables are easier to concatenate
    % to build a full data set
    pdata = struct;
    pdata.rat = config.rname;
    pdata.analysed = tnow;
    pdata.directory = pwd;
    pdata.tetrodes = tets; % list of tetrodes analysed
    pdata.sessions = nsess; % the number of recording sessions used
    pdata.data_dirs = data_dirs;   
    pdata.mapset = mapset;
    pdata.cname = config.cname;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Position data
    disp(sprintf('\t...positions')); disp(sprintf('\b '));
    pos_srate = 50; % desired sampling rate of position data (Hz)
    [pos,data_intervals,tstart] = get_pos_for_klustest('Neuralynx',data_dirs,pos_srate,mapset); % directories
    % pos = [x,y,t,v,hd,ahv] for all data
    % data_intervals = [tstart,tend] one row per recording
    % tstart = the original start time of the first recording
    pdata.pos = single( pos );
    pdata.pos_srate = pos_srate;
    total_time = sum(diff(data_intervals,[],2));
    disp(sprintf('\t\t...total recording time = %.1fs (%.1f mins)',total_time,total_time/60))
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Cluster data
    disp(sprintf('\t...clusters'))
    clus = get_clu_for_klustest('Neuralynx',config,tets);
    % clus = cell array, one cell per tetrode, cluster number of each spike
    clu_n = cellfun(@(x) numel(unique(x))-1,clus,'UniformOutput',0); % number of clusters on each tetrode
    nclus = sum([clu_n{:}]);
    disp(sprintf('\t\t...%d total clusters %s',nclus))
    pdata.clusters = clus;
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Spike data
    disp(sprintf('\t...spikes'))
    [spk,wav,spk_srate] = get_spk_for_klustest('Neuralynx',data_dirs,tets,tstart);
    spk_n = cellfun(@numel,spk,'UniformOutput',0); % number of spikes on each tetrode
    disp(sprintf('\t\t...%d spikes %s',sum([spk_n{:}])))
    pdata.spike_times = spk;
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Cluster quality
    disp(sprintf('\t...cluster quality'))
    isods = cell(1,length(tets));
    fets = cell(1,length(tets));
    quals = cell(2,length(tets));
    for tt = 1:length(tets)
        if tt==1
            disp(sprintf('\t\ttetrode: %d',tets(tt)))             
        else
            disp(sprintf('\b, %d',tets(tt))) 
        end
        fetfile = [config.cname '.fet.' num2str(tets(tt))];
        cutfile = [config.cname '_' num2str(tets(tt)) '.cut'];
        [isod,lrat,ndists,cdists,fdata] = get_iso_for_klustest('Axona',cutfile,fetfile);
        isods(1,tt) = {[isod(:) lrat(:)]};
        fets(1,tt) = {fdata};
        quals(1,tt) = {cdists};
        quals(2,tt) = {ndists};
    end
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LFP data
    disp(sprintf('\t...LFP'))
    [lfp,lfp_srate] = get_lfp_for_klustest('Neuralynx',data_dirs,tstart); % lfp = zscored amplitude, time, theta phase, theta power
    disp(sprintf('\t\t...done'))
           
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Part data
    disp(sprintf('\t...part config'))
    nparts = size(part_config,1);
    for pp = 1:nparts
        switch part_config.part_methods{pp}
            case {1} % if the method is whole session (i.e. everything recorded)
                part_config.part_times(pp) = {data_intervals(:)};
                % include all the periods where data were recorded

            case {2} % if the method is recording session
                part_config.part_times(pp) = {data_intervals(part_config.part_intervals{pp},:)};
                % include all the recording periods listed in part_intervals

            case {3} % if the method is intervals (keypresses)
                % still to add, use manageKEYS and Nlx2MatEV
                
        end
        part_config.part_duration(pp) = sum(part_config.part_times{pp}(:,2)-part_config.part_times{pp}(:,1));        
    end    
    pdata.part_config = part_config; % a copy of the part_config, this will be extended to include more data that the saved part_config though, table format
    disp(sprintf('\t\t...done'))

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Run through clusters
    disp(sprintf('Analysing clusters...'))
    sdata = table;

    % Kilosort mainly ignores tetrodes and just assigns cluster IDs across 
    % all available channels (even though we well it not to cluster across 
    % tetrodes). So, in this loop we run through every unique cluster returned
    % by kilosort and analyse its activity
    for tt = 1:length(tets) % for every tetrode
        disp(sprintf('\tTetrode %d',tets(tt)))
        clu = clus{tt}; % clusters on this tetrode
        clusters = unique(clu(clu>0)); % list on nonzero clusters
        spt = spk{tt}; % spike times for this tetrode
        wav_now = wav{tt}; % waveforms for this tetrode
        
        for cc = 1:length(clusters) % for every cluster   
            disp(sprintf('\t\tcluster %d: ',clusters(cc)))

            % The next loop focuses on dividing the data into its different parts 
            % and analysing the cluster within each of these. The results are then 
            % added to sdata_temp which will be concatenated with sdata at the end    
            for pp = 1:nparts % for every part
                sdata_temp = table; % temporary table for this cluster
                part_now = part_config.part_names{pp};
                if pp==1; disp(sprintf('\b%s ',part_now)); else; disp(sprintf('\b| %s ',part_now)); end

                % cut the position data to include only this part
                if ~isfield(pdata,part_now) || ~isfield(pdata.(part_now),'pox')
                    part_times = part_config.part_times{pp};
                    pot = pdata.pos(:,3);
                    pindax = logical(sum(pot' >= part_times(:,1) & pot' <= part_times(:,2),1));
                    ppot = pdata.pos(pindax,3); % pos time for this part            
                    ppox = pdata.pos(pindax,1); % pos x for this part
                    ppoy = pdata.pos(pindax,2); % pos y for this part
                    ppoh = pdata.pos(pindax,5); % pos HD for this part
                    ppov = pdata.pos(pindax,4); % pos HD for this part
                    ppoa = pdata.pos(pindax,6); % pos HD for this part

                    % accumulate data in pdata
                    pdata.(part_now).pox = single(ppox);
                    pdata.(part_now).poy = single(ppoy);
                    pdata.(part_now).poh = single(ppoh);
                    pdata.(part_now).pov = single(ppov);
                    pdata.(part_now).pav = single(ppoa);            
                    pdata.(part_now).pot = single(ppot);
                end

                % cut the spike data to include only this part and cluster
                sindax = logical(sum(spt' > part_times(:,1) & spt' < part_times(:,2) & (clu==clusters(cc))',1)); % spikes               
                pspt = spt(sindax); % spike time for this part
                sidx = knnsearch(ppot,pspt); % nearest neighbour in position data for every spike
                pspx = ppox(sidx); % spike x for this part
                pspy = ppoy(sidx); % spike y for this part
                psph = ppoh(sidx); % spike HD for this part
                pspv = ppov(sidx); % spike HD for this part
                pspa = ppoa(sidx); % spike HD for this part

                % accumulate data in sdata_temp
                sdata_temp.rat = {pdata.rat}; % the rat number/name
                sdata_temp.partn = uint8(pp);
    %             uci = ['r' rname '_' pdata.date '_t' num2str(tet) '_c' num2str(clu)]; % unique cell identifier - a string with [rat name, session date, electrode, cluster];
    %             sdata_temp.uci = {uci};
                sdata_temp.tetrode = uint8(tets(tt)); % the index of the cluster in Phy format (number assigned to cluster in Phy, may not start at zero, may have gaps etc)    
                sdata_temp.cluster = uint8(clusters(cc)); % the index of the cluster in Phy format (number assigned to cluster in Phy, may not start at zero, may have gaps etc)
                sdata_temp.spt_pot_index = {uint32(sidx)}; % this index can be used to get the spike values from the position data values i.e. pdata.(part_now).pox(sdata.spike_index)
                sdata_temp.spike_times = {single(pspt)}; % the actual spike times, single should be fine as spike times are sampled at 32kHz
                sdata_temp.nspikes = numel(pspx);
                sdata_temp.frate = numel(pspx) / (numel(ppot)*(1/pos_srate));           

%                 if ~any(sindax) % if there are no spikes in this part
%                     continue
%                 end

%% >>>>>>>>>> Waveform data
                nch = 4;
                mxs = NaN(nch,1);
                wav_means = NaN(nch,size(wav_now,1));
                wav_stds = NaN(nch,size(wav_now,1));  
                wav_now_clus = cell(1,4);
                for ww = 1:nch % for every channel
                    wav_now_clus{ww} = squeeze(wav_now(:,ww,sindax))'; % should be [spikes x samples]
                    wav_means(ww,:) = mean(wav_now_clus{ww},1,'omitnan');
                    wav_stds(ww,:) = std(double(wav_now_clus{ww}),[],1,'omitnan') ./ sqrt(size(wav_now_clus{ww},1));
                    mxs(ww,1) = max(wav_means(ww,:),[],'omitnan');
                end
                [~,widx] = sort(mxs,'descend'); % sort from largest > smallest waveform
                max_wav_means = wav_means(widx,:);
                max_wav_stds = wav_stds(widx,:);
                max_wav_mxs = mxs(widx);

                % use the waveform with the largest amplitude to calculate wave characteristics
                wnow = max_wav_means(1,:);
                [~,maxvalt] = max(wnow);
                wnow(1:maxvalt) = NaN;
                [minval,minvalt] = min(wnow,[],'omitnan');
                width_ms = (minvalt - maxvalt) * (1/spk_srate) * 1e03; % waveform width in ms

                % accumulate data in sdata_temp
                sdata_temp.wave_width = single(width_ms);
                sdata_temp.wave_amps = int16(mxs(1));
                sdata_temp.wave_mean = { int16(wav_means) };
                sdata_temp.wave_std = { int16(wav_stds) };
                pdata.wavtime = (( 1:size(wav_now,1) ) -8 ) * (1/spk_srate) * 1e03; % waveform samples, converted to seconds, then microseconds

%% >>>>>>>>>> Firing rate map
                % This function incorporates a number of different methods for generating firing rate maps and you can check the function for detailed descriptions of these
                % my favourite is the method proposed by Leutgeb et al. (2005) Independent Codes for Spatial and Episodic Memory in Hippocampal Neuronal Ensembles
                % https://doi.org/10.1126/science.1114037
                % although the fastest method is 'histogram' as this is just a straight binning method followed by smoothing, so it is the default
                speedlift = [];
                if isfield(pdata,part_now)
                    if isfield(pdata.(part_now),'speedlift')
                        speedlift = pdata.(part_now).speedlift;
                    end
                end
                rmset = mapset;
                pos = [ppox ppoy].*10;
                spk = [pspx pspy].*10;
                rmset.maplims = [min(pos(:,1)) min(pos(:,2)) max(pos(:,1)) max(pos(:,2))]; 
                rmset.srate = pos_srate;

                [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset,speedlift);         

                % Skaggs spatial information content (bits per second)
                % Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
                % https://onlinelibrary.wiley.com/doi/epdf/10.1002/%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K
                % Markus et al. (1994) Spatial information content and reliability of hippocampal CA1 neurons: effects of visual input
                % https://onlinelibrary.wiley.com/doi/epdf/10.1002/hipo.450040404
                pr = dwellmap ./ sum(dwellmap(:),'omitnan'); % dwell time probability
                ro = sum(ratemap(:) .* pr(:),'omitnan'); % overall firing rate
                si = sum(pr(:) .* (ratemap(:)./ro) .* log2(ratemap(:)./ro),'omitnan'); 

                % Skaggs spatial information content (bits per spike)
                % From Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
                % The information rate given by formula (1) is measured in bits per second. If it is
                % divided by the overall mean ring rate of the cell (expressed in spikes per second),
                % then a different kind of information rate is obtained, in units of bits per spike|let us
                % call it the information per spike. This is a measure of the specificity of the cell: the
                % more grandmotherish" the cell, the more information per spike. 
                sis = si ./ ro;

                % Sparsity
                % From Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
                % The sparsity measure is  an adaptation  to space of  a formula invented by Treves and Rolls (1 99 1); the adaptation measures the 
                % fraction of the environment  in which a cell  is active. Intuitively, a sparsity of, say, 0.1 means that the place field of the cell 
                % occupies 1/10 of the area the rat traverses.
                sp = (sum(pr(:).*ratemap(:),'omitnan').^2) ./ (sum(pr(:).*(ratemap(:).^2),'omitnan'));            

                % Spatial Coherence (Cacucci et al. 2007)
                % https://dx.doi.org/10.1523%2FJNEUROSCI.1704-07.2007
                % The spatial coherence for each firing rate map was computed as the mean correlation between the firing rate of each bin with the 
                % aggregate rate of the 24 nearest bins.
                meanf = ones([5 5]);
                meanf(3,3) = 0;
                U = imfilter(ratemap,meanf,0,'same','conv');
                cohe = corr(ratemap(:),U(:),'rows','pairwise','type','Pearson');            

                % accumulate data
                sdata_temp.ratemap = { single(ratemap) };
                sdata_temp.spikemap = { single(spikemap) };
                sdata_temp.spatial_info = [si sis sp cohe]; % [spatial info (b/s), spatial info (b/sp), sparsity, coherence]
                pdata.(part_now).dwellmap = dwellmap;
                pdata.(part_now).speedlift = speedlift;

%% >>>>>>>>>> Place fields
                % threshold the firing ratemap using the minimum firing rate proportion, frcut, and the minimum cutoff, minfr
                thresh_ratemap = imbinarize(ratemap,max([0, mapset.frcut*max(ratemap(:),[],'omitnan')]) );

                % run regionprops to get properties of these regions & remove fields that are too small
                datout = regionprops('table',thresh_ratemap,ratemap,'Area','Centroid','WeightedCentroid','MajorAxisLength','MinorAxisLength','Orientation','ConvexHull');
                datout.Area(:) = datout.Area(:) .* ((mapset.binsize/10)^2); % convert field area to cm2
                nindx = datout.Area(:) < mapset.arcut;
                datout(nindx,:) = [];

                % accumulate data
                sdata_temp.nfields = size(datout.Area,1);
                sdata_temp.field_data = { table2cell(datout) };

%% >>>>>>>>>> Grid score
                % create autocorrelation
                % I'm not entirely sure who to attribute this to, but autocorrelations are pretty ubiquitous now
                % this function was adapted from xPearson (adapted for speed and to make it work on n-dimensions)
                automap = ndautoCORR(ratemap,ratemap,50);

                % autocorrelation analysis
                % There are a number of different grid score approaches contained here so see the function itself for a description
                % my preferred method is the default, which is from Langston et al. (2010) Development of the Spatial Representation System in the Rat
                % In this method rings are cut from the autocorrelation at different distances from the centre and the standard rotation and correlation
                % is performed on each one. The highest of these grid scores is used, the radius of the ring is the grid spacing. A sine wave
                % is fitted to the values in the ring and this is used to estimate the positions of the fields. The grid orientation is the
                % angle from the centre to the first one of these fields, counter-clockwise. The grid field radius is taken as the radius of the centre field.
                [~,gdata] = gridSCORE(automap);            

                % accumulate data
                sdata_temp.gridmap = { single(automap) };
                sdata_temp.grid_info = single([gdata.g_score gdata.wavelength gdata.radius gdata.orientation]);

%% >>>>>>>>>> Rayleigh vector length
                % HD firing rate maps are now computed by mapHD, like the above it will accept a precomputed
                % dwell time map to speed up computation
                % rayleigh vector values etc are computed on sum normalised HD firing rate maps to remove
                % any possible effects of firing rate (it's unclear if circ_stats care about the values)
                % If you recorded HD using two LEDs this should be used here and throughout the data
                % if not, HD is estimated using the animal's displacement
                % Load a head direction dwell map if one exists, or create it if not
                hd_dwellmap = [];
                if isfield(pdata,part_now)
                    if isfield(pdata.(part_now),'hd_dwellmap')
                        hd_dwellmap = pdata.(part_now).hd_dwellmap;
                    end
                end            
                [~,hd_dwellmap,hd_spikemap,hd_ratemap,r,mx,mn,sd] = mapHD(hd_dwellmap,ppoh,psph,mapset);

                % accumulate data
                sdata_temp.hd_ratemap = { single(hd_ratemap) }; 
                sdata_temp.hd_info = single([r mx(1) mn sd]); % Rayleigh v, PFD (angle with highest firing), mean angle (norm circ mean), SD of firing (norm circ stdev)
                pdata.(part_now).hd_dwellmap = hd_dwellmap;

%% >>>>>>>>>> AHV analysis (still to add)
                poa_2cm = 2*floor(ppoa/2); % bin the ahv data in 2deg/s increments
                b_spikes = histcounts(sidx,(1:length(ppot)+1)-0.5);
                fh = fspecial('gaussian',[1 8]); % gaussian filter, 8 samples or 250ms
                insta_rate = imfilter(b_spikes,fh,'replicate','same') ./ (1/pos_srate); % smooth histogram and divide by the initial bin size

                deg_span = -360:2:360;
                ahv = NaN(length(deg_span),1);
                ahv_dwell = NaN(length(deg_span),1);            
                for dd = 1:length(deg_span)
                    ahv(dd) = mean( insta_rate(poa_2cm==deg_span(dd)),'omitnan');
                    ahv_dwell(dd) = sum(poa_2cm==deg_span(dd)) .* (1/pos_srate);
                end

                % accumulate data
                sdata_temp.ahv_curve = { single(ahv(:)) };             
                pdata.ahv_xvalues = single(deg_span);    
                pdata.(part_now).ahv_dwell = single(ahv_dwell);

%% >>>>>>>>>> Theta phase preference
                % bin the theta phase data
                ai = reshape(deg2rad(-180:5:540),[],1); % doing this means we have bins symmetrical around zero

                % we must calculate the spike phase here precisely for every spike time, rather than indexing into the position theta phase
                % this is for reasons of resolution (LFP and spikes are sampled at 48kHz positions are only sampled at 50Hz)
                % This process is fairly well established:
                % "To obtain a theta phase angle for each spike, LFPs were first bandpass filtered (fourth-order Chebyshev, r = 0.5, MATLAB 
                % filter and filtfilt routines; 6–10 Hz) before a Hilbert transform was applied to obtain the instantaneous phase angle."
                % from van der Meer and Redish (2011) Theta Phase Precession in Rat Ventral Striatum Links Place and Reward Information (https://doi.org/10.1523/JNEUROSCI.4869-10.2011)              
                pspp = interp1(lfp(:,2),lfp(:,3),pspt,'nearest');
                spp2 = [pspp(:)-2*pi; pspp(:); pspp(:)+2*pi]; % repeat phase values at +/-360 degrees
                y = histcounts(spp2,ai); % bin phase values

                % directional analyses on phase data
                mx2p = ai(y == max(y)); % preferred angle (location of max frate)

                % accumulate data
                sdata_temp.theta_phase = { single(y) }; 
                sdata_temp.theta_info = single([circ_r(pspp(:)) mx2p(1) circ_mean(pspp)]); % Rayleigh v, preferred phase (bin with max count), mean angle (average phase) 

%% >>>>>>>>>> ISI and autocorrelation analyses
                %% inter-spike interval
                [~,idata,isis] = getISIhalfwidth(pspt,100);

                % accumulate data
                sdata_temp.isi = { single(idata.adist(:,2)) }; 
                sdata_temp.isi_info = single([idata.fwhmx ]); 
                pdata.isi_xvalues = single(idata.adist(:,1));

                %% burst index
                % from Mizuseki et al. (2012) Activity Dynamics and Behavioral Correlates of CA3 and CA1 Hippocampal Pyramidal Neurons 
                % "The spike-burst index was defined as the fraction of spikes with <6 ms ISIs (Harris et al., 2001)"
                % https://doi.org/10.1002/hipo.22002  
                bindx = [isis; NaN] < 6 | [NaN; isis] < 6; % bindx is an index of all spikes sharing an isi less than 6ms
                burst_index = sum(bindx) / numel(pspt);

                % accumulate data
                sdata_temp.burst_index = single(burst_index); 

                %% Intrinsic theta    
                % from van der Meer and Redish (2011) Theta Phase Precession in Rat Ventral Striatum Links Place and Reward Information
                % "To quantify the degree and frequency of theta modulation in single cells, we used the method used by Royer et al. (2010). 
                % First we computed the autocorrelogram of the cell, in 10 ms bins from -500 to +500 ms, normalized it to the maximum value
                % between 100 and 150 ms (corresponding to theta modulation), and clipped all values above 1. Then we fit the following function [function in paper]
                % where t is the autocorrelogram time lag, and a-c, w, and t1–2 were fit using the fminsearch optimization function in MATLAB. A measure
                % of theta modulation strength, the “theta index,” was defined as a/b, which intuitively corresponds to the ratio of the sine fit relative to
                % the baseline in the autocorrelogram."
                % https://doi.org/10.1523/JNEUROSCI.4869-10.2011
                [thi,thf,thr,t500,c500,f500] = getTHETAintrinsic(pspt); % we also extract the frequency of theta and the goodness of the theta fit

                % accumulate data
                sdata_temp.autocorr_500 = { single([c500(:) f500(:)]) }; 
                sdata_temp.autocorr_500_info = single([thi thf thr]); % theta index, theta frequency, fit rsquare
                pdata.autocorr_500_xvalues = single(t500);

                %% Refractory period analyses           
                % refractory period violation analyses
                % from Navratilova et al. (2016) Grids from bands, or bands from grids? An examination of the effects of single unit contamination on grid cell firing fields.
                % although they took these analyses from pre-existing papers
                % "There are few reliable methods for estimating the contamination of a unit isolated from tetrode recordings (see DISCUSSION). The only method that does not rely 
                % on the same measures that are used for spike sorting is to check the number of spikes that occur within the refractory period of another spike..."
                % https://doi.org/10.1152/jn.00699.2015
                [nrpv,prpv,fp1,fp2,censored,t25,c25] = getRPVcount(pspt,isis,part_config.part_duration(pp)); % we extract the number of spikes in the refractory period, what proportion of the total this is etc                    

                % accumulate data
                sdata_temp.autocorr_25 = { single(c25(:)) }; 
                sdata_temp.autocorr_25_info = single([nrpv prpv]);     
                pdata.autocorr_25_xvalues = single(t25);  

%% >>>>>>>>>> Speed cell analysis
                % This analysis is taken from Kropff et al (mEC speed cell paper)
                % It is unlikely to be of great use to HPC people, although many place cells are speed modulated
                % Related to this, you can also compute the relationship between theta power and running speed
                [svals,stime,sscore,sslope,sintcpt,crve,~] = getSPEEDmod(ppot,ppov,pspt);

                % accumulate data
                sdata_temp.speed_slope = { single(interp1(svals,crve,0:1:50,'linear')) };             
                sdata_temp.speed_info = single([sscore sslope sintcpt]);             
                pdata.(part_now).speed_dwell_time = single(interp1(svals,stime,0:1:50,'linear'));

%% >>>>>>>>>> Create summary figure for this cluster in this part
    putvar(pdata,sdata_temp,wav_now_clus)
    return
                kilofig_part(pdata,sdata_temp,pp,wav_now_clus,'on');

            end % pp = 1:nparts
        end % cc = 1:nclus
    end
    
return















































