function KlusterAnalysis(rat,elecs,clus)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Main function for first pass analysis of clusters
%   KlusterAnalysis(rat,elecs,clus)
%
%%%%%%%% Inputs
%   rat = the rat's number, this can't be retrieved from anywhere so should be given as an input, this should be entered as a string
%   elecs = the electrodes to run on, given as a vector e.g. [1 2 3 4 5 6 7 8], 0 means run on 1-8 (def = all recorded tetrodes)
%   clus = the clusters to run on, set this to 0 to run on all clusters (def = all)
%
%%%%%%%% Comments
%   23/02/16 YTcluanalysisROD majorly revised
%   24/02/16 added functionality to find session details from .set file
%   30/03/16 fixed some problems with .set file data because I updated it to handle multiple sessions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get starting parameters
function_started = datestr(now);
disp(sprintf('Starting at: %s',function_started));  
disp(sprintf('\t...extracting session data: merge_set.mat'));
load('merge_set.mat'); % Extract setdata information from merged .set file

% Get a list of tetrodes that were used
if ~exist('elecs','var') || isempty(elecs)
    elecs = setdata.active_tetrodes;
elseif elecs == 0
    elecs = [1 2 3 4 5 6 7 8];
end % if isempty(elecs) || ~exist('elecs','var')

% Define which clusters to run on
if ~exist('clus','var') || isempty(clus)
    clus = 0;
    clus_ref = clus;
else
    clus_ref = clus;
end % if isempty(elecs) || ~exist('elecs','var')

% Define the name of the session to analyse - KlustaKwiker should output any analysis as 'merge' files
if isfield(setdata,'out_name')
    filename = setdata.out_name;
else
    filename = 'merge';
end % if isfield(sdata.out_name)

% Get some important sampling information
pixel_ratio = setdata.pixels_per_metre;  
dt_position = setdata.sample_rate; % position sampling rate (ms)

%% Display starting messages                                              
disp(sprintf('\t...running on: %s',pwd));
disp(sprintf('\t...processing: %s',filename));

% electrodes and clusters
disp(sprintf('\t...analysing tetrodes: %s',mat2str(elecs)));
if clus ~= 0
    disp(sprintf('\t...analysing clusters: %s',mat2str(clus)));
else
    disp(sprintf('\t...analysing all clusters'));
end % if clus ~= 0

% Digital inputs
if setdata.digital_inputs == 1 && exist([filename '.key'],'file')
    is_key_divide = 1;
    disp(sprintf('\t...using digital inputs: %s',[filename '.key']));
else
    is_key_divide = 0;
    disp(sprintf('\t...no digital inputs found'));
end % if is_key_divide == 1 && exist(sdata.filenames.key)

ka_start = tic; % Start timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional Parameters
fformat = '-dpng';
fquality = '-r300';
figvis = 'off';

bin_size = 2.5; % (cm), for calculating the rate map.
sigma = 1.5; % sigma (gaussian standard deviation) to be used for rate and position map smoothing 
min_dwell = 0.00001;                                                                                                                          % total number of seconds that rat has to be in a bin for it to count

smoo = 1.5; % smoothing (sigma) of maps
fsiz = [5 5]; % smoothing window of maps (length x height), must be odd numbers
mdwell = 0.05; % minimum dwell time for dwell time map

overwrite = 0; % set to 1 if you want KlusterData file to be overwritten by whatever is produced in this analysis (the function will move the last KlusterData to Old_KlusterData, but this will be overwritten the next time it happens)

% output xml configuration file (Klusters will not work without this)
write_config(filename,50,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Position Data
[filenames,posx,posy,hd,post,timeseg,seg_size] = read_mypos(filename); % Position data should already be analysed by KlusterKwiker
shift_array = []; 
if isempty(hd) || mean(hd) == 0
    disp(sprintf('Actual HD not found...'));
    if exist('MatFiles\KlusterData.mat');
        disp(sprintf('\t...loading estimated HD'));
        load('MatFiles\KlusterData.mat');
        struct_date = ['s_' setdata.s_1.converted_date];
        struct_rat = ['r_' rat];
        hd = KlusterData.(struct_rat).(struct_date).head_direction;
    else
        disp(sprintf('\t...estimating HD'));
        [hd] = estimateHD(posx,posy,dt_position,1,5);
    end % if exist('MatFiles\KlusterData.mat');
end % if isempty(hd) || mean(hd) == 0
disp(sprintf('\t...done'));

try                                                                                                                                     % The next line will fail if a the .goal file is not consistent with the data - so I added a nice explanatory error. Roddy.
    time_index = ones(1, length(seg_size));
    for ii=2:length(time_index)
            time_index(ii) = time_index(ii-1) + seg_size(ii-1);
    end % for ii=2:length(time_index)
catch errgoal2
    disp(sprintf('\tERROR: KlusterAnalysis has encountered an error. Check that your .goal file is correct...'));
    return
end % try/catch loop

session_names = {};
session_ind = 1;
bucket_session_names = {};
bucket_ind = 1;
for ii=1:length(filenames)
    fname = filenames{ii};
    ind = regexp(upper(fname),'BK','Once');
    if(isempty(ind))
        session_names{session_ind} = fname;
        session_ind = session_ind+1;
    else
        bucket_session_names{bucket_ind} = fname;
        bucket_ind = bucket_ind+1;
    end % if(isempty(ind))
end % for ii=1:length(filenames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup sep_config
% sep_config = {} means separating all merged files into separate analyses output
% sep_config = 0 means producing one big merged analysis output
% sep_config = n means producing one output for the selected trajectory
if ~exist([filename '.key'],'file') || ~exist([filename '.goal'],'file')
	sep_config = 0;
    % configuring default sep_config conditions
	if isempty(sep_config)
		for ii=1:length(time_index)
		        sep_config{ii} = ii;
		end % for ii=1:length(time_index)
        out_names = filenames;
	elseif(~iscell(sep_config))
		if sep_config == 0
		        sep_config = {};
		        sep_config{1} = 1:length(time_index);
		        out_names{1} = filename;
		else
		        tmp = sep_config;
		        sep_config = {};
		        sep_config{1} = tmp;
		        out_names{1} = filenames{tmp};
		end % if(sep_config ==0)
	end % if(isempty(sep_config))
else
	% Cluanalysis now uses the .goal file
    sep_config = {};
	goals = dlmread(goals_filename);						
	goal_array = [];
	for ii=1:length(goals(:,1))
		goal_array = [goal_array goals(ii,goals(ii,:)>0)];
	end % for ii=1:length(goals(:,1))

	sep_config{1} = find(goal_array == 1);								                                        % Find all 1's in .goal file (runs to goal 1)  
	sep_config{2} = find(goal_array == 2); 								                                        % Find all 1's in .goal file (runs to goal 1)
	sep_config{3} = find(goal_array == 3);								                                        % Find all 2's in .goal file (runs to goal 2)
	sep_config{4} = find(goal_array == 4);								                                        % Find all 3's in .goal file (runs to goal 3)
	sep_config{5} = [sep_config{1} sep_config{2} sep_config{3} sep_config{4}];								                                        % Find all 4's in .goal file (runs to goal 4)
	out_names={'Cylinder 1' 'Parallel' 'Radial' 'Cylinder 2' 'Entire Session'}; 

	if length(sep_config) ~= length(out_names)
		disp(sprintf('WARNING: sep_configs and outnames do not match, exiting'));
		return
	end % if length(sep_config) ~= length(outnames)
end % if screening == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Spike Data
if overwrite
    if ~exist('MatFiles\KlusterData.mat')
        movefile('MatFiles\KlusterData.mat','MatFiles\Old_KlusterData.mat')
        delete('MatFiles\KlusterData.mat')
    end % if ~exist('MatFiles\KlusterData.mat') || overwrite;
end % if overwrite

% as the recording gains should not change between merged sessions, choose the first file to represent
disp('Reading gains...');
setfn = [filenames{1} '.SET'];
base_volts = gain_info(setfn);
disp(sprintf('\t...done')); 

orig_posx = posx;
orig_posy = posy;
orig_hd = hd;
posx = round(posx);
posy = round(posy);
hd = round(hd);

disp('Starting analysis...');
disp(sprintf('\t...session lasted %.f minutes (%.fs)',setdata.total_session_length/60,setdata.total_session_length));

[~,~,~] = mkdir('MatFiles'); % KlusterData will be stored here
[~,~,~] = mkdir('Figures'); % Figures will be stored here
for e = 1:length(elecs)
    eno = elecs(e);
    fetfile = [filename '.fet.' num2str(eno)]; %out_names={'Goal 1' 'Goal 2' 'Goal 3' 'Goal 4' 'Entire Maze' 'Square Box' }; 
    clufile = [filename '.clu.' num2str(eno)];
    spkfile = [filename '.spk.' num2str(eno)];

    % check number of clusters
    [nclu,labels] = getclusters2(clufile);
    disp(sprintf('Analysing tetrode %d...', eno));

    if (nclu == 1)
            disp(sprintf('\t...done'));
            disp(sprintf('\t...electrode %d: found 0 clusters, skipped...', eno));                   
    else
        % load stuff
        fid = fopen(spkfile);
        if (fid == -1)
            disp(sprintf('\tERROR: File not found: %s',spkfile));
        else
            waves = fread(fid,Inf,'int16=>int8');
            waves = reshape(waves, 4, 50, []);
        end %  if (fid == -1)
        fclose(fid);

        fid = fopen(fetfile);
        if (fid == -1)
            disp(sprintf('\tERROR: File not found: %s',fetfile));
            continue;
        else
            nfeatures = sscanf(fgetl(fid),'%d');
            features = fscanf(fid,'%d',[nfeatures Inf]);
            fclose(fid);
        end %  if (fid == -1)

        features = int32(features);
        fet_per_ted = (nfeatures -4) / 4;
        num_spikes = length(features(1,:));
        disp(sprintf('\t...done'));                                                                                                             % display 'done' - indented in the command window
        disp(sprintf('\t...electrode %d: found %d clusters, %d spikes total.',eno,nclu-1,length(labels)));                 % Print in a nice green

        % calculate cluster quality
        % only uses energy and PCA for calculating L-ratio
        [isods,lratios,cd2s,nd2s] = cluquality(features, nclu, labels, [1 3 1+fet_per_ted 3+fet_per_ted 1+fet_per_ted*2 3+fet_per_ted*2 1+fet_per_ted*3 3+fet_per_ted*3]);
        color_pars = choose_color(nclu);

        if (clus_ref(1) == 0)
            clus = [];
            for cc =1:nclu
                clus = [clus cc];
            end % for cc =1:nclu
        end % if (clus_ref(1) == 0)

        % open figure 2 - the figure containing cluster energy plots                
        fig_cspace = figure('Visible',figvis);
        rowa1 = 4; % Number of rows for subplots
        colu1 = 3; % Number of columns for subplots
        sa_spac1 = 0.03; % Spacing between plots
        sa_padd1 = 0.03;
        sa_marg1 = 0.01;
        fsize1 = 4;
        axlwid = 1;
        for cc = 1:length(clus)
            clu = clus(cc);
            disp(sprintf('\t\t...cluster %d',clu-1));

            fets = features(:,find(labels == clu));
            waves_c = waves(:,:,find(labels == clu));
            hds = fets(nfeatures-3,:);
            posxs = fets(nfeatures-2,:);
            posys = fets(nfeatures-1,:);
            spikest = fets(nfeatures,:);
            spikest = double(spikest) ./ 100000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Cluster Data
            color_par = color_pars(clu,:);
            fet_vec = [1 6; 1 11; 1 16; 6 11; 6 16; 11 16];
            lab_vec = {{'energy1' 'energy2'} {'energy1' 'energy3'} {'energy1' 'energy4'} {'energy2' 'energy3'} {'energy2' 'energy4'} {'energy3' 'energy4'}};
            for ff = 1:length(fet_vec(:,1))
                fet1 = fet_vec(ff,1);
                fet2 = fet_vec(ff,2);
                lablz = lab_vec{ff};
                labl1 = lablz{1};
                labl2 = lablz{2};
                
                set(0, 'CurrentFigure', fig_cspace); % set current figure to fig_two without changing its visibility setting
                subaxis(rowa1,colu1,ff,'Spacing',sa_spac1,'Padding',sa_padd1,'Margin',sa_marg1);
                hold on
                h1 = plot(fets(fet1,:),fets(fet2,:),'.','MarkerSize',2);
                xlabel(labl1);
                ylabel(labl2);
                set(h1,'Color',color_par);
                
                title([labl1 ' vs ' labl2]);
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize1);
                box off
                axis square
                
                subaxis(rowa1,colu1,ff+6,'Spacing',sa_spac1,'Padding',sa_padd1,'Margin',sa_marg1);
                hold on
                h1 = plot(fets(fet1,:),fets(fet2,:),'.','MarkerSize',1);
                xlabel(labl1);
                ylabel(labl2);
                set(h1,'Color',[0.5 0.5 0.5]);
                title([labl1 ' vs ' labl2]);
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize1);
                box off 
                hold on
                
                rings = 5; % number of rings to show, each ring represents 1 standard deviation
                COMx = nanmean(nanmean(fets(fet1,:)));
                COMy = nanmean(nanmean(fets(fet2,:)));
                dists = fastDist([COMx COMy],[fets(fet1,:)' fets(fet2,:)']);
                std_dist = nanstd(dists);
                
                cents = repmat([COMx COMy],rings,1);
                radis = repmat(std_dist,rings,1) .* (1:rings)';

                viscircles(cents,radis,'EdgeColor',color_par,'LineWidth',0.5);
                
                title([labl1 ' vs ' labl2]);
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize1);
                box off
                axis square
            end % for ff = 1:length(fet_vec(:,1))
            
% Ignore noise cluster (which is conventionally assigned to cluster 1)
            if (clu == 1) 
                continue;
            end % if (clu == 1)

% Calculate root mean squared values for channel noise - for calculating signal to noise ratio later
% SNR as seen in Stratton, Cheung, Wiles, Kiyatkin, Sah, et al. (2012) Action Potential Waveform Variability...
            waveforms_n = waves(:,:,find(labels == 1));
            waveforms_n = double(waveforms_n);
% Convert the waveform value to the actual microvolts
            waveforms_n(1,:,:) = waveforms_n(1,:,:)*base_volts(4*(e-1) + 1);
            waveforms_n(2,:,:) = waveforms_n(2,:,:)*base_volts(4*(e-1) + 2);
            waveforms_n(3,:,:) = waveforms_n(3,:,:)*base_volts(4*(e-1) + 3);
            waveforms_n(4,:,:) = waveforms_n(4,:,:)*base_volts(4*(e-1) + 4);
% Calculate root mean squares (RMS)
            rms_noise1 = mean(mean(rms(squeeze(waveforms_n(1,:,:))))); % rms channel 1
            rms_noise2 = mean(mean(rms(squeeze(waveforms_n(2,:,:))))); % rms channel 2
            rms_noise3 = mean(mean(rms(squeeze(waveforms_n(3,:,:))))); % rms channel 3
            rms_noise4 = mean(mean(rms(squeeze(waveforms_n(4,:,:))))); % rms channel 4
            rmss = [rms_noise1 rms_noise2 rms_noise3 rms_noise4]; % accumulate root mean squares for later
                     
% sep_config = the goals - in the order set up at the start
% sepi will = which goal we are looking at, e.g 1 for goal one                       
% start separating the merged sessions
            for sepi=1:length(sep_config)
                try                                                                                                     % Cluanalysis will fail here if the out_names are not set up properly, so I added an error message. Roddy.
                    id = [out_names{sepi} '-elec' num2str(eno) 'clu' num2str(clu)];
                    sep_array = sep_config{sepi};
                    sep_array = sort(sep_array);
                    nspikes = 0;
                    total_time = 0;
                    goodspikes = 0;
                catch errsepi
                    disp(sprintf('\tERROR: Cluanalysis has encountered an error. Check that out_names at start of script are correct (i.e. right number of out_names)'));
                    return
                end % try/catch loop

% holders for the position information
                sep_pos_index = ones(1,length(sep_array));
                posx_now = [];
                posy_now = [];
                hd_now = [];
                post_now = [];

% holders for the spikes information
                sep_spk_index = ones(1,length(sep_array));
                cspikes = [];
                waves_ct = [];
                posxst = [];
                posyst = [];
                hdst = [];

% sep_array = which trajectories correspond to the goal we are looking at
% ari will = the trajectory number e.g. 45 for the 45th trial/trajectory
                for ari=1:length(sep_array)
                    segno = sep_array(ari);
                    start_ti = time_index(segno); % if error: check your .goal file is correct
                    stop_ti = 0;

                    tstart = post(time_index(segno));
                    if (segno == length(timeseg))
                        tstop = post(end) + 5;
                        stop_ti = length(post);
                    else
                        tstop = post(time_index(segno+1));
                        stop_ti = time_index(segno+1) - 1;
                    end % if (segno == length(timeseg))

% retrieve the session related information
                    posx_now_part = posx(start_ti:stop_ti);
                    posy_now_part = posy(start_ti:stop_ti);
                    hd_now_part = hd(start_ti:stop_ti);
                    post_now_part = post(start_ti:stop_ti);

                    posx_now = [posx_now; posx_now_part];
                    posy_now = [posy_now; posy_now_part];
                    post_now = [post_now; post_now_part];
                    hd_now = [hd_now; hd_now_part];
                    if (ari < length(sep_array))
                        sep_pos_index(ari+1) = sep_pos_index(ari) + length(post_now_part);
                    end % if (ari < length(sep_array))

% get spike times
                    cspikes_part = spikest(spikest >= tstart);
                    waves_ct_part = waves_c(:,:,spikest>=tstart);
                    hdst_part = hds(spikest >= tstart);
                    posxst_part = posxs(spikest >= tstart);
                    posyst_part = posys(spikest >= tstart);

                    cspikest = spikest(spikest >= tstart);

                    cspikes_part = cspikes_part(cspikest < tstop);
                    waves_ct_part = waves_ct_part(:,:,cspikest < tstop);
                    hdst_part = hdst_part(cspikest < tstop);
                    posxst_part = posxst_part(cspikest < tstop);
                    posyst_part = posyst_part(cspikest < tstop);

                    nspikes = nspikes + length(cspikes_part);
                    total_time = total_time + (tstop-tstart);
                    avg_trial_time = total_time/length(sep_array);

                    spkhd = hdst_part(find(~isnan(hdst_part)));
                    goodspikes = goodspikes + length(spkhd);

                    cspikes = [cspikes; cspikes_part'];
                    posxst = [posxst; posxst_part'];
                    posyst = [posyst; posyst_part'];
                    hdst = [hdst; hdst_part'];
                    waves_ct = cat(3, waves_ct, waves_ct_part);
                    if (ari < length(sep_array))
                            sep_spk_index(ari+1) = sep_spk_index(ari) + length(cspikes_part);
                    end % if (ari < length(sep_array))
                end % for ari=1:length(sep_array)

                indx = find(posxst > nanmax(posx_now) | posxst < nanmin(posx_now));
                posxst(indx,:) = [];
                posyst(indx,:) = [];
                indx = find(posyst > nanmax(posy_now) | posyst < nanmin(posy_now));
                posxst(indx,:) = [];
                posyst(indx,:) = [];

                % Determine spike frequency and skip if necessary
                if(length(posxst) ==0)
                    continue;
                end  %if(length(posxst) ==0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings for figure 1 plots
                fig_cluster = figure('Visible',figvis);
                rowa = 3; % Number of rows for subplots
                colu = 4; % Number of columns for subplots
                sa_spac = 0.015; % Spacing between plots
                sa_padd = 0.02;
                sa_marg = 0.03;
                fsize2 = 4;
                
                % tile locations for fig_cluster
                spike_tile = 1; % spikes and trajectory
                dwell_tile = 2; % dwell time map
                frate_tile = 3; % firing rate map
                grids_tile = 4; % grid autocorrelogram tile
                headd_tile = 5; % head direction tile
                theta_tile = 6; % theta locking tile
                autor_tile = 7; % autocorrelogram (refractory period) tile
                autot_tile = 8; % autocorrelogram (theta) tile
                wavef_tile = 9; % waveform tile
                spikt_tile = 10; % spikes over time tile
                mehal_tile = 11; % spikes over time tile
                thfit_tile = 12; % Theta fitted line tile
                
                %id = ['E' num2str(eno) 'C' num2str(clu) ' - YTcluanalysis - ' out_names{sepi}];
                %annotation('textbox', [0.05, 1.0, 1.0, 0], 'string', ([sprintf('%s: %d/%d spks (%.2f Hz) Spike Peak %.1f Spike Width %d',id,goodspikes,nspikes,nspikes/total_time, max_amp, spk_width)]),'FontSize',10);                    % [x y w h]                                                         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: Spikes vs Time
                subaxis(rowa,colu,spikt_tile,'Spacing',sa_spac,'Padding',sa_padd,'Margin',sa_marg,'Holdaxis');
                tstart = min(post_now);
                tstop = max(post_now);
                tstop = tstop + 5; % padding at the end for stray spikes
                timeh = hist(cspikes,ceil(tstop - tstart));
                bar(floor(tstart):(floor(tstart) + ceil(tstop -tstart) -1), timeh,0.9,'k');
                ylabel('Frate (Hz)');
                if (timeh <= 0)
                        timeh = 1;
                end %  if (timeh <= 0)

                axis([tstart tstop 0 max(timeh)*1.1]);
                xlabel('Time (s)');
                box off
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: Place field (firing rate map) plot
                subaxis(rowa,colu,frate_tile,'Spacing',sa_spac,'Padding',sa_padd,'Margin',sa_marg,'Holdaxis');
                ind = ~isnan(posx_now);
                posx_ = posx_now(ind);
                posy_ = posy_now(ind);
                post_ = post_now(ind);

                % get rid of data points where both posx and posy = 0
                sum_pos = posx_ + posy_;
                ind = find(sum_pos>0);
                posx_ = posx_(ind);
                posy_ = posy_(ind);
                post_ = post_(ind);
                skaggs = NaN;
                spars = NaN;
                cohe = NaN;
                frmap = [];

                % prepare data for mapping
                posx_ = double(posx_);
                posy_ = double(posy_);
                orig_hd = double(orig_hd);
                posxst = double(posxst);
                posyst = double(posyst);

                % Generate map
                miny1 = min(-posy_);
                minx1 = min(posx_);
                
                pox = posx_ - minx1 + 1;
                poy = -posy_ - miny1 + 1;
                pot = post_;
                spx = posxst - minx1 + 1;
                spy = -posyst - miny1 + 1;

                maxxnow = ceil(nanmax(pox));
                maxynow = ceil(nanmax(poy));

                [mapout1,~] = getMAP(pox,poy,maxxnow,maxynow,pixel_ratio,bin_size); % dwell time map
                [mapout2,~] = getMAP(spx,spy,maxxnow,maxynow,pixel_ratio,bin_size); % spike map

                mapout1 = imgaussfilt(mapout1,smoo,'FilterSize',fsiz);
                mapout2 = imgaussfilt(mapout2,smoo,'FilterSize',fsiz);

                mapout1 = mapout1 .* dt_position/1000;
                mapout1(mapout1 < mdwell) = NaN;
                frmap = mapout2 ./ mapout1;

                indx = isnan(frmap);
                frmap(indx) = 0;
                frmap = imgaussfilt(frmap,1.5,'FilterSize',[5 5]);
                frmap(indx) = NaN;

                posmap = mapout1;
                spikmap = mapout2;
                
                skaggs = skaggs_info2(frmap,posmap);
                spars = sparsity(frmap, posmap);
                cohe = spatial_coherence(frmap, posmap);
                
                %% Calculate shuffled spatial information content
                datin = struct;
                datin.pox = pox;
                datin.poy = poy;
                datin.pot = pot;
                datin.spx = spx;
                datin.spy = spy;
                datin.sot = cspikes;
                datin.maxxnow = maxxnow;
                datin.maxynow = maxynow;
                datin.pixel_ratio = pixel_ratio;
                datin.bin_size = bin_size;
                datin.smoo = smoo;
                datin.fsiz = fsiz;
                datin.mdwell = mdwell;
                datin.dt_position = dt_position;
                datin.shuffn = 100;
%                 [datout] = skaggsSHUFFLED(datin);
% 
%                 
%                 dbstop
                
                
                
                
                
                
                
                % Plot map
                im = imagesc(frmap);
                set(im,'alphadata',~isnan(frmap))
                colormap(jet(256));
                axis xy
                axis off
                daspect([1 1 1])
                title(sprintf('SI %.2f b/s   SP %.2f%%   Cohe %.2f',skaggs, (spars*100), cohe));
                caxis([0 nanmax(nanmax(frmap))]);
                
                % Add colorbar and reduce its width while retaining the position of both
                x1 = get(gca,'position');
                c = colorbar;
                set(gca,'position',x1)
                x = get(c,'Position');
                x(3) = 0.01;
                set(c,'Position',x)
                
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: Spatial autocorrelation plot (for grid cells)
                subaxis(rowa,colu,grids_tile,'Spacing',sa_spac,'Padding',sa_padd,'Margin',sa_marg,'Holdaxis');
                AutoCorrM = GridAutoCorr(frmap);
                [grid_score,grid_spacing,field_size,grid_orientation,grid_ellipticity] = GridAnalysis(AutoCorrM,bin_size);

                % Plot map
                AutoCorrM = AutoCorrM ./ nanmax(nanmax(AutoCorrM));
                im = imagesc(AutoCorrM);
                set(im,'alphadata',~isnan(AutoCorrM))
                colormap(jet(256));
                axis xy
                axis off
                daspect([1 1 1])
                if isnan(grid_score)
                    title(sprintf('Too few fields to calculate'));
                else
                    title(sprintf('G: %.1f Spacing: %.1f Size: %.1f Ori: %.1f Elli: %.1f',grid_score,grid_spacing,field_size,grid_orientation,grid_ellipticity));
                end % if isnan(grid_score)
                caxis([-1 1]);
                
                % Add colorbar and reduce its width while retaining the position of both
                x1 = get(gca,'position');
                c = colorbar;
                set(gca,'position',x1)
                x = get(c,'Position');
                x(3) = 0.01;
                set(c,'Position',x)
                
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: Chronographic (dwell time) plot
                subaxis(rowa,colu,dwell_tile,'Spacing',sa_spac,'Padding',sa_padd,'Margin',sa_marg,'Holdaxis');
                % Plot map
                im = imagesc(posmap);
                set(im,'alphadata',~isnan(posmap))
                colormap(jet(256));
                axis xy
                axis off
                daspect([1 1 1])
                title(sprintf('Time: %.2f s',total_time));
                caxis([0 nanmax(nanmax(posmap))])
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize2);

                % Add colorbar and reduce its width while retaining the position of both
                x1 = get(gca,'position');
                c = colorbar;
                set(gca,'position',x1);
                x = get(c,'Position');
                x(3) = 0.01;
                set(c,'Position',x);
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: Position plot
                subaxis(rowa,colu,spike_tile,'Spacing',sa_spac,'Padding',sa_padd,'Margin',sa_marg,'Holdaxis');
                hold on;
                % Plot position data points
                for ari=1:length(sep_pos_index)
                    if (ari == length(sep_pos_index))
                        plot(pox(sep_pos_index(ari):end), poy(sep_pos_index(ari):end),'Color',[0.5 0.5 0.5]);
                    else
                        plot(pox(sep_pos_index(ari):(sep_pos_index(ari+1)-1)), poy(sep_pos_index(ari):(sep_pos_index(ari+1)-1)),'Color',[0.5 0.5 0.5]);
                    end % if (ari == length(sep_pos_index))
                end % for ari=1:length(sep_pos_index)

                % Plot spike data points
                for ari=1:length(sep_spk_index)
                    if (ari == length(sep_spk_index))
                        plot(spx(sep_spk_index(ari):end), spy(sep_spk_index(ari):end), 'r.', 'MarkerSize',10);
                    else
                        plot(spx(sep_spk_index(ari):(sep_spk_index(ari+1)-1)), spy(sep_spk_index(ari):(sep_spk_index(ari+1)-1)), 'r.', 'MarkerSize',10);
                    end % if (ari == length(sep_spk_index))
                end % for ari=1:length(sep_spk_index)

                axis([min(min(pox)) max(max(pox)) min(min(poy)) max(max(poy))]);
                set(gca, 'DataAspectRatio', [1 1 1]);
                title(sprintf('Spikes: %.f',nspikes));
                axis off;
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize2);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: Alternative Head Direction (HD) plot - line graph rather than polar plot
                bins = (0:5:360)';
                [hdfr,ahb,shb,spike_hd,rayleigh_v,pfd_mean,pfd_std,pfd_max,max_fr] = analyseHD(hd,bins,dt_position,cspikes);

                subaxis(rowa,colu,headd_tile,'Spacing',sa_spac,'Padding',sa_padd,'Margin',sa_marg,'Holdaxis');
                [ax,p1,p2] = plotyy(bins,hdfr,bins,ahb,'plot','plot');
                xlabel(ax(1),'Head Direction'); % label x-axis
                ylabel(ax(1),'Firing Rate (Hz)'); % label left y-axis
                ylabel(ax(2),'Occupancy (s)'); % label right y-axis
                p1.LineWidth = 2;
                p1.Color = 'k';

                set(ax(1),'LineWidth',axlwid);
                p2.LineWidth = 1;
                p2.Color = [0.5 0.5 0.5];
                set(ax(2),'LineWidth',axlwid);
                set(ax,{'ycolor'},{'k'; [0.5 0.5 0.5]});
                box off
                set(ax(1),'xlim',[0 360]);
                set(ax(2),'xlim',[0 360]);
                title(sprintf('r = %.2f, mu = %.2f, max = %.2f',rayleigh_v,pfd_mean,pfd_max));
                set(ax(1),'LineWidth',axlwid,'layer','top','FontSize',fsize2);
                set(ax(2),'LineWidth',axlwid,'layer','top','FontSize',fsize2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: Mehalanobis Distance        
                subaxis(rowa,colu,mehal_tile,'Spacing',sa_spac,'Padding',sa_padd,'Margin',sa_marg,'Holdaxis');
                nd2 = nd2s{clu};
                cd2 = cd2s{clu};
                max_d = max([max(nd2);max(cd2)]);
                max_d = round(max_d);

                % to cope with missing channel data
                if (isnan(max_d)) 
                        max_d = 0;
                end % if (isnan(max_d))

                hv = 0:1:max_d;
                isod = isods(clu);
                lratio = lratios(clu);

                [cdist,tbins1] = hist((cd2),hv);
                [ndist,tbins2] = hist((nd2),hv);
                plot(tbins1,cdist,'b','LineWidth',2);
                hold on;
                plot(tbins2,ndist,'k','LineWidth',2);
                set(gca,'XScale','log');
                set(gca, 'Xlim', [0,max_d + 100]);
                title(sprintf('Iso-D = %.2f, Lratio = %.2f',isod,lratio));

                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize2);
                box off
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: 200ms autocorrelation - Theta plot                                  
                subaxis(rowa,colu,autot_tile,'Spacing',sa_spac,'Padding',sa_padd,'Margin',sa_marg,'Holdaxis');
                [corrdata,tms] = spk_acorr2(cspikes,1,500); % get spike autocorrelogram
                [thetaR,thetaP,thetaIndx,thetaPowr,thetaLin] = getTHETAfit([tms' corrdata]); % fit a decomposing sine wave to the spike autocorrelogram
                 
                area(tms,corrdata,'FaceColor',[0.5 0.5 0.5]);
                hold on
                plot(tms,thetaLin,'r','linewidth',2);
                hold on
                title(sprintf('R = %.2f, p = %d, Tindx = %.2f, Tpowr = %d',thetaR,thetaP,thetaIndx,thetaPowr));
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize2);
                box off
                v1 = axis;
                axis([tms(1) tms(end) 0 v1(4)]);
                xlabel('Time Lag (ms)');
                ylabel('Spike Probability');
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: 20ms autocorrelation - Refractory period plot
                subaxis(rowa,colu,autor_tile,'Spacing',sa_spac,'Padding',sa_padd,'Margin',sa_marg,'Holdaxis');
                [corrdata,tms] = spk_acorr2(cspikes,1,50);
                bar(tms,corrdata,0.9,'k');
                %title([id ': autocorrelation']);
                v = axis;
                axis([tms(1) tms(end) v(3) v(4)]);
                xlabel('Time Lag (ms)');
                ylabel('Spike Probability');
                box off
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize2);

                %% Calculate refractory contamination according to Navratilova and McNaughton (2016) Grids from bands, or bands from grids? An examination of the effects of single unit contamination on grid cell firing fields
                half_spike = corrdata(tms >= 0); % take only the positive side of spike autocorrelogram
                half_time = tms(tms >= 0); % take only the positive side of spike autocorrelogram

                rlength = 2; % length of refractory period in ms
                lockout = 0.75; % lockout period of recording system before another spike can be detected in ms

                pviol = sum((half_spike ~= 0) & (half_time' <= rlength)); % get the number of refractory period violations
                F = nspikes/total_time;
                refrac_contam = (1-(sqrt(1-((4 * pviol)/(F * (2*(rlength-lockout)))))))/2;
                title(sprintf('Rviol = %.2f, Contamination = %.2f',pviol,refrac_contam));

                hold on
                plot([-rlength; -rlength],[0 v(4)],'r','LineWidth',1)
                plot([rlength; rlength],[0 v(4)],'r','LineWidth',1)          
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the four waveform channels
                if length(posxst) > 0
                    waveforms = waves_ct;
                    waveforms = double(waveforms);

                    waveforms(1,:,:) = waveforms(1,:,:)*base_volts(4*(e-1) + 1);						% convert the waveform value to the actual microvolts
                    waveforms(2,:,:) = waveforms(2,:,:)*base_volts(4*(e-1) + 2);
                    waveforms(3,:,:) = waveforms(3,:,:)*base_volts(4*(e-1) + 3);
                    waveforms(4,:,:) = waveforms(4,:,:)*base_volts(4*(e-1) + 4);

                    mean1 = mean(waveforms,3);

                    std1 = std(waveforms,0,3);
                    all_mean = reshape(mean1, [],1);
                    maxwav = max(all_mean)*2;
                    minwav = min(all_mean)*2;

                    rand_index = rand(50,1);
                    rand_index = ceil(rand_index*length(waveforms(1,1,:)));
                    waveforms = waveforms(:,:,rand_index);

                    channel_status_all = ones(1,4);
                    channel_s2n_all = NaN(1,4);
                    channel_peak_all = NaN(1,4);
                    channel_width_all = NaN(1,4);
                    for i = 1:4
                        rms_signal = mean(rms(squeeze(waveforms(i,:,:)))); % calcluate max peak amplitude for signal to noise ratio calculation
                        rms_noise = rmss(i); % get noise rms from accumulated vector
                        sig2noise = rms_signal/rms_noise; % calculate signal to noise ratio for this channel

                        max_signal = max(max(squeeze(mean1(i,:))));
                        if max_signal == inf
                            max_signal = max(max(squeeze(mean1(i,:))<inf));
                        end % if max_signal == inf

                        ch = squeeze(mean1(i,:));
                        [maxval,maxt] = max(ch);
                        [postminval,postmint] = min(ch(maxt:end));
                        postmint = postmint + maxt - 1;
                        width = postmint - maxt;
                        width = width * (1000/50); 

                        if mean(waveforms(i,:,:)) ==0 % if the waveforms = o the channel is dead
                            channel_status_all(i) = 0;
                        else
                            channel_s2n_all(i) = sig2noise;
                            channel_peak_all(i) = max_signal;
                            channel_width_all(i) = width;
                        end % if mean(waveforms(1,:,:)) ==0	                                        
                    end % for i = 1:4
                end % if length(spx) > 0
     
                wave_peaks = channel_peak_all;
                wave_widths = channel_width_all;
                indx = find(wave_peaks == max(max(wave_peaks)));
                max_width = wave_widths(indx);    
                    
                subaxis(rowa,colu,wavef_tile,'Spacing',sa_spac,'Padding',sa_padd,'Margin',sa_marg,'Holdaxis');
                wavtime = -200:20:780;
                plot(wavtime,squeeze(waveforms(indx,:,:)),'k','LineWidth',0.5);
                hold on;
                plot(wavtime,squeeze(mean1(indx,:))+squeeze(std1(indx,:)),'r--','LineWidth',1);
                plot(wavtime,squeeze(mean1(indx,:))-squeeze(std1(indx,:)),'r--','LineWidth',1);
                plot(wavtime,squeeze(mean1(indx,:)),'r','LineWidth',1);
                axis([-200 780 minwav maxwav]);                                                                       
                text(350,maxwav-maxwav/12,sprintf('S=%.2f,P=%.2f,W=%.2f',channel_s2n_all(indx),channel_peak_all(indx),channel_width_all(indx)),'BackgroundColor','w','FontSize',fsize2);
                box off
                set(gca,'LineWidth',axlwid,'layer','top','FontSize',fsize2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell Typing         
                frate2 = nspikes/total_time;
                wave_peaks = channel_peak_all;
                wave_widths = channel_width_all;
                indx = find(wave_peaks == max(max(wave_peaks)));
                max_width = wave_widths(indx);

                [ctype] = getCTYPE(frate2,max_width,skaggs,grid_score,rayleigh_v);       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect session data
                if ~exist('MatFiles\KlusterData.mat')
                    KlusterData = struct;
                else
                    load('MatFiles\KlusterData.mat');
                end % if ~exist('MatFiles\KlusterData.mat') || overwrite;

                s_date = setdata.s_1.converted_date;
                cUID = [rat '_' s_date '_' num2str(eno) '_' num2str(clu)]; % create unique cell id from rat no, date, electrode and cluster
                struct_rat = ['r_' rat];
                struct_date = ['s_' s_date];
                struct_eno = ['e_' num2str(eno)];
                struct_clu = ['c_' num2str(clu)];
                %% Session data
                KlusterData.(struct_rat).(struct_date).positions = [posx_ posy_];
                KlusterData.(struct_rat).(struct_date).head_direction = hd;
                %% General data
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).Unique_identity = cUID;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).KlusterKwiked = setdata.last_klustering;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).KlusterAnalysed = function_started;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).rat = rat;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).date = s_date;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).electrode = eno;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).cluster = clu;
                %% Cluster data
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).lratio = lratios(clu);
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).isolation_distance = isods(clu);
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).cell_type = ctype;
                %% Spike train data
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).refractory_period_contamination = refrac_contam;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).theta_correlation = thetaR;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).theta_p_value = thetaP;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).theta_index = thetaIndx;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).theta_power = thetaPowr;
                %% Spike Data
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).start_time = tstart;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).stop_time = tstop;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).total_time = total_time;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).total_spikes = nspikes;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).average_firing = nspikes/total_time;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).spike_hist = timeh;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).spike_times = cspikes;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).spike_autocorrelation = [corrdata tms'];
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).spikes = [posxst posyst];
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).spike_hd = spike_hd;
                %% Map data
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).ratemap = frmap;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).dwellmap = posmap;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).spikemap = spikmap;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).spatial_information = skaggs;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).sparsity = spars;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).coherence = cohe;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).max_firing = nanmax(nanmax(frmap));
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).mean_firing = nanmean(nanmean(frmap));
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).min_firing = nanmin(nanmin(frmap));
                %% Grid data
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).grid_autocorr = AutoCorrM;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).grid_score = grid_score;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).grid_spacing = grid_spacing;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).grid_size = field_size;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).grid_orientation = grid_orientation;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).grid_ellipticity = grid_ellipticity;
                %% Head direction data
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).rayleigh_vector = rayleigh_v;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).pfd_mean = pfd_mean;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).pfd_max = pfd_max;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).pfd_deviation = pfd_std;
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).max_hd_frate = max_fr;
                %% Waveform data
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).channel_status = channel_status_all;% waveforms active or dead (1/0)
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).channel_signal_to_noise = channel_s2n_all;% waveform signal to noise ratios
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).channel_peak_amplitude = channel_peak_all;% waveform peaks (uV)
                KlusterData.(struct_rat).(struct_date).(struct_eno).(struct_clu).(out_names{sepi}).channel_width_waveform = channel_width_all;% waveform widths (us)

                %% Save data file
                save('MatFiles\KlusterData.mat','KlusterData');
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save figures
                % Save the summary figure for this session/cluster (actually the second figure opened)
                set(0, 'CurrentFigure', fig_cluster);
                antext = [cUID ' ' ctype];
                annotation('textbox', [0.05, 1.0, 1.0, 0], 'string',antext,'FontSize',5);                    % [x y w h]                                                         
                %saveas(fig_cluster,['Figures\' cUID out_names{sepi}],fformat);
                print(fig_cluster,['Figures\' cUID out_names{sepi}],fformat,fquality)
                close(fig_cluster);
            end % for sepi=1:length(sep_config)
        end % for cc = 1:length(clus)

        % Save cluster space figure - acutally the first figure opened
        set(0, 'CurrentFigure', fig_cspace);
        annotation('textbox', [0.05, 1.0, 1.0, 0], 'string', (sprintf('%s: %s %s %s',id,num2str(nclu-1),'clusters',num2str(num_spikes))),'FontSize',5);                    % [x y w h]                                                         
        %saveas(fig_cspace,['Figures\Electrode ' num2str(eno) ' ClusterSpace' out_names{sepi}],fformat);
        print(fig_cspace,['Figures\Electrode ' num2str(eno) ' ClusterSpace' out_names{sepi}],fformat,fquality)
        close(fig_cspace);
    end % if (nclu == 1)
end % for e = 1:length(elecs)

%% Save data file
save('MatFiles\KlusterData.mat','KlusterData');

%% Finish up
toc1 = toc(ka_start)/60;
disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
disp(sprintf('KlusterAnalaysis has finished. It took %0.3f seconds or %0.3f minutes', toc(ka_start),toc1)) % Stop counting time and display results


















































