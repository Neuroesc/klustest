function [dat,lfp_srate] = get_lfp_for_klustest(dataformat,data_dirs,tstart,ds)
% get_lfp_for_klustest load local field potential for klustest
% Data loading function for klustest, given a data format type, and some
% settings, finds .eeg files (Tint format) or .ncs files (Neuralynx format) and
% loads the local field potential data. Can also downsample to a lower sampling
% rate, which is useful when loading raw files.
%
% USAGE
%
% [dat,lfp_srate] = get_lfp_for_klustest(dataformat,data_dirs,tstart,ds)
%
% INPUT
%
% 'dataformat' - String: 'kwikcut','neuralynx','kwiktint','phy' (experimental)
%
% 'data_dirs' - Data directories, output from get_tets_for_klustest
%
% 'tstart' - Start time of the recordings, given by get_pos_for_klustest
%
% 'ds' - Optional, new samppling rate to resample LFP to
%
% OUTPUT
%
% 'dat' - Nx4, where N = number of LFP samples, columns = LFP amplitude in
%       volts, time point of each sample, theta phase, theta power
%
% 'lfp_srate' - Sample rate of LFP, if ds was given, this will equal ds,
%               otherwise gives the native sample rate
%
% NOTES
% 1. Phy data format is experimental
%
% 2. Function is not really intended to be used without klustest
%
% 3. Theta phase and power estimated using Hilbert transform of the LFP filtered
% between 6 and 12 Hz.
%
% 4. Resampling achieved using Matlab's resample, method = 'pchip'
% 
% SEE ALSO kwiktint klustest get_lfp_volts Nlx2MatCSC

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
    if ~exist('ds','var') || isempty(ds) || all(isnan(ds))
       ds = 0;
    end
    if ~exist('tstart','var') || isempty(tstart) || all(isnan(tstart))
       tstart = 0;
    end
    
    switch lower(dataformat)
        case {'kwiktint'}
            lfpt = [];
            lfp = [];
            % lfp_channel = 1;
            total_duration = 0;
    
            for ff = 1:length(data_dirs) % for every data directory  
                [~,b,~] = fileparts(data_dirs{ff});          
                fname_eeg = [b '.eeg']; % the file name 
                % fname_egf = [b '\' b '.egf']; % the file name 
    
                if ds
                    [lfp_now,t_now,lfp_srate,lfp_h,~] = get_lfp_volts(fname_eeg,'ds',ds);
                else
                    [lfp_now,t_now,lfp_srate,lfp_h,~] = get_lfp_volts(fname_eeg);
                end      
                t_now = t_now+total_duration;
                total_duration = total_duration + lfp_h.duration;
    
                % concatenate
                lfpt = [lfpt(:); t_now(:)]; % convert microseconds to seconds and concatenate
                lfp = [lfp(:); lfp_now(:)]; % concatenate amplitudes
            end
            lfpt = lfpt - tstart;
    
            % filter LFP to get theta and calculate theta phase
            [b,a] = butter(2,[6 12]/(lfp_srate/2)); % Generate 4th order butterworth filter coefficients for theta band [6 12] Hz
            lftheta = filtfilt(b,a,lfp); % Apply filter to data using zero-phase filtering  
            h = hilbert(lftheta); % hilbert transform of the theta signal
            theta_phase = mod(angle(h),2*pi); % get the phase of theta
            theta_power = abs(h); % get the amplitude of theta (theta power) this is not normalized
            
            dat = [lfp(:) lfpt(:) theta_phase(:) theta_power(:)];
    
        case {'neuralynx'}
            lfpt = [];
            lfp = [];
            lfp_channel = 1;
            for ff = 1:length(data_dirs) % for every data directory     
                lfp_fname = [data_dirs{ff} '\CSC' num2str(lfp_channel) '.ncs'];        
                [lfpt_n, ch_n, Fs_n, ~, lfp_n, h_n] = Nlx2MatCSC(lfp_fname, [1 1 1 1 1],  1, 1, 1);  
                lfp_srate = Fs_n(1);
    
                lfp_now = lfp_n(:);
                % lfp_n is 512 x N where N is the number of records and the number
                % of timestamps in lfpt_n
                % basically, the data are split into columns of 512 samples and
                % only timestamps for the fist sample in each column is given
                t_now = (0:length(lfp_now)-1) .* (1/lfp_srate) .* 1e06 + lfpt_n(1); 
                t_now = t_now ./ 1e06; % convert microseconds to seconds
                % t is a vector same length as LPF, spaced according to sampling rate, converted to microseconds, shifted to start at beginning of session
    
                if ds
                    [lfp_now,t_now] = resample(lfp_now,t_now,ds,'pchip'); % use interpolation and an anti-aliasing filter to resample the signal at a uniform sample rate
                    lfp_srate = ds; % the new sampling rate
                end      
                
                % concatenate time
                lfpt = [lfpt(:); t_now(:)]; % convert microseconds to seconds and concatenate
    
                % convert to microvolts and concatenate LFP
                adbv_idx = contains(h_n,{'ADBitVolts'}); % find the ADBitVolts entry
                adbv = str2double( h_n{adbv_idx}(13:end) ); % extract the number part of the text and convert to double     
                lfp_now = lfp_now .* adbv .* 10^6; % convert spike AD to volts, then to microvolts
                lfp = [lfp(:); lfp_now(:)]; % concatenate amplitudes
            end
            lfpt = lfpt - tstart;
    
            % filter LFP to get theta and calculate theta phase
            [b,a] = butter(2,[6 12]/(lfp_srate/2)); % Generate 4th order butterworth filter coefficients for theta band [6 12] Hz
            lftheta = filtfilt(b,a,lfp); % Apply filter to data using zero-phase filtering  
            h = hilbert(lftheta); % hilbert transform of the theta signal
            theta_phase = mod(angle(h),2*pi); % get the phase of theta
            theta_power = abs(h); % get the amplitude of theta (theta power) this is not normalized
            
            dat = [lfp(:) lfpt(:) theta_phase(:) theta_power(:)];
            
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
    
        case {'ElePhy'}
            lfpt = [];
            lfp = [];
            dat = [];
            lfp_channel = 1;
            lfp_srate = 30e3;
            for ff = 1:length(data_dirs) % for every data directory           
    %             lfp_fname = [data_dirs{ff} '\' data_dirs{ff} '.LFP\' data_dirs{ff} '.timestamps.dat'];
    %             dat = readTrodesExtractedDataFile(lfp_fname);
            end   
    end




































