function [dat,lfp_srate] = get_lfp_for_klustest(dataformat,data_dirs,tstart,ds)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  short desc.
% long description
%
% USAGE:
%       [out] = template(in) process with default settings
% 
%       [out] = template(in,optional1) process using optional argument 1
% 
%       [out] = template(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%       of the process
% 
%       Parameters include:
% 
%       'param1'          -   (default = X) Scalar value, parameter to do something
% 
%       'param2'          -   (default = X) Scalar value, parameter to do something
% 
% INPUT:
%       in    - input as a vector
%
% OUTPUT:
%       out   - output as a vector
%
% EXAMPLES:
%       % run function using default values
%       out = template(in,varargin)
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
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
if ~exist('ds','var') || isempty(ds) || all(isnan(ds))
   ds = 0;
end
if ~exist('tstart','var') || isempty(tstart) || all(isnan(tstart))
   tstart = 0;
end

switch dataformat
    case {'Neuralynx'}
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




































