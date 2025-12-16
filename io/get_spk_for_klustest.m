function [spk,wav,spk_srate,wavtime] = get_spk_for_klustest(dataformat,data_dirs,tets,tstart,clus,spike_window)
% get_spk_for_klustest load spike waveforms and spike times
% Data loading function for klustest, given a data format type, and some
% settings, loads waveforms and spike times
%
% USAGE
%
% [spk,wav,spk_srate,wavtime] = get_spk_for_klustest(dataformat,data_dirs,tets,tstart,clus,spike_window)
%
% INPUT
%
% 'dataformat' - String: 'kwikcut','neuralynx','kwiktint','phy' (experimental)
%
% 'data_dirs' - Data directories, output from get_tets_for_klustest
%
% 'tets' - Vector, electrodes we want to analyse
%
% 'tstart' - Start time of the recordings, given by get_pos_for_klustest
%
% 'clus' - Vector, Nx1 where N = number of spikes, gives the cluster assigned to
%       each spike, given by get_clu_for_klustest
%
% 'spike_window' - [-0.2 0.8], the time window over which to load and plot waveforms (only used when loading 
%       Phy data, but used for all plotting)
%
% OUTPUT
%
% 'spk' - Cell array, spike times in seconds for all spikes, each cell is an
%       electrode
%
% 'wav' - Cell array, waveforms for all spikes, each cell is an electrode
%       (samples x channels x spikes)
%
% 'spk_srate' - sample rate of the waveforms in Hz
%
% 'wavtime' - time in us of each waveform sample relative to the peak
%
% NOTES
% 1. Phy data format is experimental
%
% 2. Function is not really intended to be used without klustest
% 
% SEE ALSO kwiktint klustest get_dacq_headers

% HISTORY
%
% version 1.0.0, Release 24/08/16 Initial release, decided to create specific loading functions for klustest
% version 1.0.1, Release 04/04/18 Commenting and formatting for klustest update
% version 1.0.2, Release 29/11/25 Updates for GitHub release
% version 1.0.3, Release 16/12/25 updated filenames for cross platform flexibility
% version 1.0.4, Release 16/12/25 improved comments, not ideal but detailed enough for now
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2025 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    switch dataformat
        case {'kwiktint'}
            spk_srate = 48e3;
            spk = cell(1,size(tets,1));
            wav = cell(1,size(tets,1));
    
            for tt = 1:size(tets,1) % for every available electrode
                all_spt = [];
                all_wav = [];
                total_duration = 0;
    
                for ff = 1:length(data_dirs) % for every data directory 
                    [~,b,~] = fileparts(data_dirs{ff});          
                    fname = [b '.' num2str(tets(tt,1))]; % the file name                
    
                    [spt,ch1,ch2,ch3,ch4,set_h] = get_spikes_volts(fname);
                    h = get_dacq_headers(fname); % sometimes the setfile doesn't contain the duration
                    spt = spt+total_duration;
                    total_duration = total_duration + h.duration;
                    all_spt = [all_spt; spt(:)];
    
                    % waveforms
                    spikes = permute(cat(3,ch1,ch2,ch3,ch4),[2 3 1]); % samps x channels x spikes
                    all_wav = cat(3,all_wav,spikes);
                end
                spk{tt} = all_spt - tstart;
                wav{tt} = all_wav;
                % wavtime = (( 1:size(spikes,1) ) -8 ) * (1/spk_srate) * 1e03;
                wavtime = linspace(0,1000,50)-200;
                % From the dacqUSB manual (pp 40)
                % In single unit mode, the system records only a 1 ms threshold-triggered "window" 
                % of samples for each spike (data stored for 200 s pre- and 800 s post-threshold), 
                % along with a timestamp indicating when the spike occurred relative to the beginning 
                % of the trial.
    
                disp(sprintf('\t\t...%d spikes %s',numel(all_spt)))
            end              
    
        case {'kwikcut'}    
            spk_srate = 32e3;
            spk = cell(1,size(tets,1));
            wav = cell(1,size(tets,1));
    
            for tt = 1:size(tets,1) % for every available electrode
                all_spt = [];
                all_wav = [];
                for ff = 1:length(data_dirs) % for every data directory 
                    electrode_type = tets(tt,2);
                    [a,b,c] = fileparts(data_dirs{ff});                
                    if electrode_type==0 % tetrode
                        fname = fullfile(b,['TT' num2str(tets(tt,1)) '.ntt']); % the file name                
                    elseif electrode_type==1 % stereotrode
                        fname = fullfile(b,['ST' num2str(tets(tt,1)) '.nst']); % the file name                                
                    end
                    [spt,spikes,h] = Nlx2MatSpike(fname,[1 0 0 0 1],1,1,1);                 
                    spt = ( spt / 1e6 ) - tstart; % convert from microseconds to seconds, then rezero to start of first recording
                    all_spt = [all_spt; spt(:)];
                    
                    adbv_idx = contains(h,{'ADBitVolts'}); % find the ADBitVolts entry
                    adbv_txt = h{adbv_idx}(13:end); % extract the number part of the text
                    adbv = str2double(strsplit(adbv_txt)); % split and convert to double
                    spikes = spikes .* adbv .* 10^6; % convert spike AD to volts, then to microvolts
                    all_wav = cat(3,all_wav,spikes);
                end
                spk{tt} = all_spt;
                wav{tt} = all_wav;
                wavtime = (( 1:size(spikes,1) ) -8 ) * (1/spk_srate) * 1e03;
    
                disp(sprintf('\t\t...%d spikes %s',numel(all_spt)))
            end              
    
        case {'neuralynx'}
            spk_srate = 32e3;
            spk = cell(1,size(tets,1));
            wav = cell(1,size(tets,1));
            for tt = 1:size(tets,1) % for every available electrode
                all_spt = [];
                all_wav = [];
                for ff = 1:length(data_dirs) % for every data directory     
                    [spt, ~, ~, ~, spikes, h] = Nlx2MatSpike( fullfile(data_dirs{ff},['TT' num2str(tets(tt,1)) '.ntt']),ones(1,5),1,1,1); 
                    spt = ( spt / 1e6 ) - tstart; % convert from microseconds to seconds, then rezero to start of first recording
                    all_spt = [all_spt; spt(:)];
                    
                    adbv_idx = contains(h,{'ADBitVolts'}); % find the ADBitVolts entry
                    adbv_txt = h{adbv_idx}(13:end); % extract the number part of the text
                    adbv = str2double(strsplit(adbv_txt)); % split and convert to double
                    spikes = spikes .* adbv .* 10^6; % convert spike AD to volts, then to microvolts
                    all_wav = cat(3,all_wav,spikes);
                end
                spk{tt} = all_spt;
                wav{tt} = all_wav;
                wavtime = (( 1:size(spikes,1) ) -8 ) * (1/spk_srate) * 1e03;
    
                disp(sprintf('\t\t...%d spikes %s',numel(all_spt)))
            end
            
        case {'phy'}
            spk_srate = 30e3;
            spk = cell(1,max(tets));
            wav = cell(1,max(tets));
            for tt = 1:length(tets) % for every available electrode
                all_spt = [];
                all_wav = [];
                for ff = 1:length(data_dirs) % for every data directory         
                    spike_times_fname = fullfile(data_dirs{ff},[data_dirs{ff} '.mountainsort'],['output_T' num2str(tets(tt))],'phy_MS','spike_times.npy');
                    spk_idx = readNPY(spike_times_fname); % all spike times for all tetrodes
                    spt = single(spk_idx) .* (1/spk_srate); % convert from samples from start to seconds from start
                    all_spt = [all_spt; spt(:)];
                    
                    gwfparams.dataDir = fullfile(data_dirs{ff},[data_dirs{ff} '.mountainsort'],['output_T' num2str(tets(tt))],'phy_MS');    % KiloSort/Phy output folder
                    gwfparams.fileName = 'recording.dat'; % .dat file containing the raw waveforms
                    gwfparams.dataType = 'int16'; % Data type of .dat file (this should be BP filtered)
                    gwfparams.nCh = 4; % Number of channels that were streamed to disk in .dat file
                    gwfparams.wfWin = round(spike_window*1e-03*spk_srate); % Number of samples before and after spiketime to include in waveform
                    clu = clus{tets(tt)};
                    gwfparams.spikeTimes = spk_idx(clu>0); % Vector of cluster spike times (in samples) same length as .spikeClusters
                    gwfparams.spikeClusters = clu(clu>0); % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
                    
                    % Load .dat and KiloSort/Phy output
                    fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);           
                    filenamestruct = dir(fileName);
                    dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
                    nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
                    wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
                    mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'}); % memmap the data file
                    nch = 4; % number of recording channels
                    waves = NaN(wfNSamples,nch,numel(spk_idx));
                    idx = find(clu>0);
                    spk_idx2 = spk_idx(idx);
                    
                    for ww = 1:numel(spk_idx2) % for every data spike
                        wav_now = mmf.Data.x(1:nch,spk_idx2(ww)+gwfparams.wfWin(1):spk_idx2(ww)+gwfparams.wfWin(end)); % extract the spikes for all 4 channels
                        waves(:,:,idx(ww)) = wav_now'; % put them in waves
                    end                
                    all_wav = cat(3,all_wav,-waves);     
    
                end
                spk{tets(tt)} = all_spt;
                wav{tets(tt)} = all_wav;
                wavtime = (gwfparams.wfWin(1) : 1 : gwfparams.wfWin(2)) * (1/spk_srate) * 1e03;
    
                disp(sprintf('\t\t...%d spikes %s',numel(all_spt)))            
            end
            
    end

% Note from Neuralynx:
% Neuralynx stores all data in our record files as A/D values. This means that any program you use that reads these files, 
% will read the data in as A/D values by default, NOT volts. In order to convert these A/D values to volts, you need to 
% view the header of the file in question and find the “-ADBitVolts” entry. This value will give you the number of volts 
% for each bit increment in the A/D value. To convert to volts, simply multiply the A/D value by the ADBitVolts value.
% 
% Voltage(V)=ADValue*ADBitVolts.
% 
% Some Neuralynx programs (such as Neuraview) have an option to display data values as either A/D values or volts. See 
% the Neuraview manual for more information. The input range of the Digital Lynx is +/- 132mV





















































