function [lfpv,t,Fs,lfp_h,set_h] = get_lfp_volts(fname,varargin)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% get_lfp_volts  function to import Axona LPF data and convert to microvolts
% This function acts as a wrapper for get_lfp, whilst also converting the LFP data 
% to microvolts using the method outlined at the end of this file (supplied by Jim)
% This function can also downsample LFP data if required (useful for huge EGF files)
%
% USAGE:
%       [lfpv,t,Fs,h] = get_lfp_volts(fname) process with default settings
% 
%       [lfpv,t,Fs,h] = get_lfp_volts(___,Name,Value,...) process with Name-Value pairs used to control aspects of the process
% 
%       Parameters include:
% 
%       'ds'   -   (default = 0) Scalar value, downsample data to this sample rate (e.g. 100 = downsample to 100Hz)
% 
% INPUT:
%       fname    -  String, lfp filename, can be full file path, .eeg or .egf
%
% OUTPUT:
%       lfpv   - Vector, Nx1 where N = number of samples, local field potential converted to microvolts
%
%       t   - Vector, Nx1 where N = number of samples, time of each sample
%
%       Fs   - Scalar, sample rate of lfp (will equal ds if provided as an input)
%
%       lfp_h   - Structure, headers for LFP file
%
%       set_h   - Structure, headers for .set file
%
% EXAMPLES:
%       % run function, downsample to 100 Hz
%       [lfpv,t,Fs,h] = get_lfp_volts('tint_data.eeg','ds',100)
%
% See also: get_lfp FUNCTION3

% HISTORY:
% version 1.0.0, Release 27/04/18 Initial release
% version 1.1.0, Release 07/08/18 updated and cleaned
% version 2.0.0, Release 07/08/18 incorporated getLFP3 here so getLFPV can do everything
% version 2.0.1, Release 07/08/18 fixed bug in the way .EEG files are converted
% version 3.0.0, Release 19/10/21 renamed get_lfp_volts, simplified, changed to use get_dacq_headers and get_lfp
% version 3.0.1, Release 19/10/21 removed need for .set file name
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
    addRequired(p,'fname',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    addParameter(p,'ds',0,@(x) isnumeric(x) && isscalar(x));   
    parse(p,fname,varargin{:});
    config = p.Results;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % load the LFP data from the file
    [lfp,t,Fs,lfp_h] = get_lfp(config.fname,'ds',config.ds);

    % Convert LFP - work out which eeg channel this is
    [a,b,c] = fileparts(config.fname);
    ech = regexp(c,'\d*','Match');
    if isempty(ech)
        ech = 1;
    else
        ech = str2double(ech{1});
    end

    % get the channel this EEG channel corresponds to (i.e. EEG channel 2 might be recording channel 34 etc)
    set_file_name = [a '\' b '.set'];    
    if ~exist(set_file_name,'file')
        error(sprintf('\t...file not found: %s',set_file_name))
    end
    set_h = get_dacq_headers(set_file_name);    
    rch = set_h.(['EEG_ch_' num2str(ech)]);

    % get the gain for this channel and the system ADC
    gain = set_h.(['gain_ch_' num2str(rch-1)]);
    adc = set_h.ADC_fullscale_mv;

    % convert LFP to microvolts
    if strcmp(c,'.eeg')
        lfpv = lfp ./ 128 .* adc .* 1000 ./ gain; % these are stored at 8-bit resolution, so the values range from -128 to 127
    else
        lfpv = lfp ./ 32768 .* adc .* 1000 ./ gain; % these are stored at 16-bit resolution, so the values range from -32768 to 32767
    end

% Raw text from Jim outlining the process above:
%
%     The values in the files are in bits.  That is, they have not been adjusted
%     for the gain factor.  To convert those numbers to microvolts, you need to
%     do three things:
% 
%     1. Look up the gain for the channel in the .set file for the trial.  The
%     relevant values are called gain_ch_X, where X goes from 0 to 127 for
%     channels 1 to 128.  So, if you want to convert the samples for channel 1
%     to volts, look up gain_ch_0.
% 
%     2. Next, look up the value of ADC_fullscale_mv in the .set file.  It is
%     almost always 1500 for dacqUSB systems.
% 
%     3. The data values for .egf files range from -32768 to 32767 (16 bits).
% 
%     To convert from sample value to microvolts, the factor is:
% 
%     value_in_uV = sample_value / 32768 * ADC_fullscale_mv * 1000 / gain
% 
%     where "sample_value" is the signed-integer value you read from the file,
%     "ADC_fullscale_mv" is the value you found at step 2 above, and "gain" is
%     the relevant gain factor you found at step 1.
% 
%     For spike (.N) or .eeg data, these are stored at 8-bit resolution, so the
%     values range from -128 to 127, and the conversion factor is:
% 
%     value_in_uV = sample_value / 128 * ADC_fullscale_mv * 1000 / gain














