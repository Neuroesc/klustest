 function [lfp,t,Fs,h,d] = get_lfp(fname,varargin)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% get_lfp  load LFP data from Axona format files
% loads LFP data from Axona files, can downsample data if necessary
% both .eeg and .egf files are acceptable
%
% USAGE:
%       [lfp,t,Fs] = get_lfp(fname) process with default settings
% 
%       [lfp,t,Fs] = get_lfp(___,Name,Value,...) process with Name-Value pairs used to control aspects of the process
% 
%       Parameters include:
% 
%       'ds'   -   (default = 0) Scalar value, downsample data to this sample rate (e.g. 100 = downsample to 100Hz)
% 
% INPUT:
%       fname    -  String, filename, can be full file path
%
% OUTPUT:
%       lfp   - Vector, Nx1 where N = number of samples, local field potential
%
%       t   - Vector, Nx1 where N = number of samples, time of each sample
%
%       Fs   - Scalar, sample rate of lfp (will equal ds if provided as an input)
%
%       h   - Structure, headers for LFP file
%
%       d   - Scalar, line where data starts
%
% EXAMPLES:
%       % run function using default values
%       [lfp,t,Fs] = get_lfp('tint_data.eeg','ds',150)
%
% See also: get_dacq_headers resample

% HISTORY:
% version 1.0.0, Release 26/02/16 initial release of getEEG and getEGF
% version 2.0.0, Release 11/01/17 created as a faster and more compact alternative to getEEG and getEGF
% version 2.1.0, Release 01/08/17 took advantage of some shortcuts in a function Jim sent, added downsampling
% version 2.2.0, Release 07/08/18 updated and cleaned
% version 3.0.0, Release 19/10/21 renamed get_lfp, changed to use get_dacq_headers, added comments
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    p = inputParser;
    addRequired(p,'fname',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    addParameter(p,'ds',0,@(x) isnumeric(x) && isscalar(x));   
    parse(p,fname,varargin{:});
    config = p.Results;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % get header info and data starting line
    [h,d] = get_dacq_headers(fname);
    Fs = h.sample_rate;
    nsamp = h.num_EEG_samples;

    % move forward to 'data start' position
    fid = fopen(fname,'r','ieee-be'); % open file    
    fseek(fid,d,0);

    % read in the LFP data
    [~,~,c] = fileparts(fname);
    if strcmp(c,'.eeg')
        lfp = fread(fid,nsamp,'int8');
    else
        lfp = fread(fid,nsamp,'int16');
    end    
    t = (0:length(lfp)-1)'/Fs;
    fclose(fid);

    % downsample data if required
    if config.ds
        [lfp,t] = resample(lfp,t,config.ds,'pchip'); % use interpolation and an anti-aliasing filter to resample the signal at a uniform sample rate
        Fs = config.ds; % the new sampling rate
    end
    lfp = double(lfp(:));
    t = double(t(:));
    Fs = double(Fs);














