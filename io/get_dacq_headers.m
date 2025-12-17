function [h,d] = get_dacq_headers(fname)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% get_dacq_headers  retrieves headers from Axona file formats
% This function reads the headers contained at the top of different Axona files
% They are added to a structure where the structure field names are the header
% names (i.e. 'trial_time' or 'sample_rate') and the values are those given after
% the header. The values are tested to see if they are a) a string, b) a number or
% c) a number followed by a units string like 'Hz'. These are left as strings (a) 
% or converted to double floats (b) with the units removed (c).
%
% USAGE:
%       [h,d] = get_dacq_headers(fname) process with default settings
% 
% INPUT:
%       fname - String, filename, can be full file path
%
% OUTPUT:
%       h  - structure, contains headers read from file
%       e.g. fields may include but are not limited to:
%         h.trial_time = string corresponding to trial onset
%         h.experimenter = text from dacqUSB 'experimenter' text box
%         h.comments = text from dacqUSB 'comments' text box
%         h.duration = duration of recording in s
%         h.num_colours = number of colours tracked
%         h.min_x = minimum x value of camera window, max_x, min_y, max_y also exist
%         h.window_min_x = minimum x value of tracking window, window_max_x, window_min_y, window_max_y also exist
%         h.timebase = frequency of position sampling
%         h.pixels_per_metre = pixels per metre value taken from dacqUSB 'pixel ratio' text box
%         h.num_pos_samples = the number of position data samples recorded, should be duration*timebase
%
%       d - line number where headers end and data actually starts
%
% EXAMPLES:
%       % run function using default values
%       h = get_dacq_headers('tint_data.eeg');
%
% See also: fopen

% HISTORY:
% version 1.0.0, Release 21/07/17 created to replace completely stupid headerCheck and key_value functions, to allow variable pixel ratio etc
% version 1.1.0, Release 21/07/17 finished .pos and .set conditions
% version 1.2.0, Release 24/07/17 added LED colours and bearings from .set file
% version 1.3.0, Release 24/07/17 added functionality for all file types, realised that they can all be contained in the same loop
% version 1.4.0, Release 24/07/17 regexp has some buggy effects, swicthed to strfind, streamlined conversion of numbers to digits
% version 2.0.0, Release 19/10/21 renamed get_dacq_headers, uses a comparison list and general file reading to collect headers
% version 2.0.1, Release 19/10/21 added comments, tested on different files, working
% version 2.0.1, Release 19/10/21 removed need for comparison lists, function is fully automatic
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    h = struct;
    d = 0;
    fid = fopen(fname,'r','ieee-be'); % open file
    if fid < 0
        return
    end    

    line = fgets(fid); % look at the next line of the file (end line = -1)
    h = struct;
    ind = 1;
    while line ~= -1 % run through every line until the end, theoretically
        if length(line)>9 & contains('data_start',line(1:10)) 
            d = ind;
            break
        else
            % break the line into the header name, which is the first string (i.e. 'trial_date')
            % and the second part which is the header info (i.e. 'Tuesday, 18 Oct 2021')
            [t,r] = strtok(line);
            
            % We test if the header info is completely numeric, if it is we convert to a double float
            % if not we will test if it is a number followed by a space and some units (i.e. '50 Hz')
            % If this is true, we will convert the first number to a double float and drop the units
            r = strtrim( r(2:end) ); % remove whitespaces and return characters at each end of the header info
            if isnan( str2double(r) ) % if the header info contains characters
                c = strsplit(r,' '); % split the info into parts based on whitespaces (i.e. '50 Hz' becomes '50' and 'Hz')
                if isnan( str2double(c{1}) ) % if the info still contains text it can't simply be units, so keep it as a string
                    h.( t ) = r;
                else % if we are just left with numbers it is likely there was just a units string on the end, so convert to double
                    h.( t ) = str2double(c{1});                        
                end
            else % if the header info is all numbers
                h.( t ) = str2double(r);
            end

            % get next line
            line = fgets(fid);  
        end     
        ind = ind+1;
    end    
    fclose(fid);


































