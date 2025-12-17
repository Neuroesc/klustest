function [ts,ch1,ch2,ch3,ch4,set_h] = get_spikes_volts(fname,varargin)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% get_spikes_volts  function to import Axona spike data and convert to microvolts
% This function acts as a wrapper for get_spikes, whilst also converting the data 
% to microvolts using the method outlined at the end of this file (supplied by Jim)
%
% USAGE:
%       [ts,ch1,ch2,ch3,ch4,set_h] = get_spikes_volts(fname) process with default settings
% 
%       [ts,ch1,ch2,ch3,ch4,set_h] = get_spikes_volts(fname___,Name,Value,...) process with Name-Value pairs used to control aspects 
%       of the process
% 
%       Parameters include:
% 
%       'bfast'  -   (default = 0) Scalar value, if 1 the function will only output the spike times, ts, the other outputs will be set to NaN
% 
% INPUT:
%       fname    -  String, lfp filename, can be full file path, .eeg or .egf
%
% OUTPUT:
%       ts       - 1x50 vector, the spike times (in seconds)
%
%       ch1-ch4  - Mx50 single matrix, the spike waveforms in microvolts, each row is a spike each column is a sample
%
%       set_h    - Structure, headers for .set file
%
% EXAMPLES:
%       % run function using default values
%       [ts,ch1,ch2,ch3,ch4,set_h] = get_spikes_volts('tint_data.1','bfast',1)
%
% See also: get_spikes get_dacq_headers

% HISTORY:
% version 1.0.0, Release 02/04/19 Initial release
% version 2.0.0, Release 19/10/21 get_spikes replaces getspikes, get_dacq_headers replaces getDACQHEADERS
% version 2.0.1, Release 19/10/21 updated for non-planar grid project
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
    addParameter(p,'bfast',0,@(x) isnumeric(x) && isscalar(x));   
    parse(p,fname,varargin{:});
    config = p.Results;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % load the spike data from the file
    [ts,ch1,ch2,ch3,ch4] = get_spikes(config.fname,config.bfast);
    if config.bfast % if bfast, we don't want waveforms, ch1-ch4 will be NaN and conversion to uV is not necessary
        return
    end

    % Convert spikes - work out which tetrode this is
    [a,b,c] = fileparts(fname);
    if isempty(a)
        a = pwd;
    end
    tet = str2double(regexp(c,'\d+','match'));

    % get the data for this session
    set_file_name = [a '\' b '.set'];
    if ~exist(set_file_name,'file')
        error(sprintf('\t...file not found: %s',set_file_name))
    end
    set_h = get_dacq_headers(set_file_name);   

    % work out the channels corresponding to this tetrode
    chindx = reshape((0:127)',4,[])';
    chindx = chindx(tet,:);

    % get their gains and convert waveforms to microvolts
    adc = set_h.ADC_fullscale_mv;
    chs = {ch1 ch2 ch3 ch4};
    for gg = 1:length(chindx)
        gain = set_h.(['gain_ch_' num2str(chindx(gg))]);
        chs{gg} = single(chs{gg}) ./ 128 .* adc .* 1000 ./ gain;
    end
    ch1 = chs{1};
    ch2 = chs{2};
    ch3 = chs{3};
    ch4 = chs{4};

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














