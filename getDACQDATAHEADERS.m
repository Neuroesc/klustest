function hdata = getDACQDATAHEADERS(snames,ftype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function can read the header info from any dacq file type (i.e. spike file, .pos, .eeg, .set)
%   and outputs all the data as a structure, where rows correspond to sessions, field names equal the 
%   variable name in the dacq file and the value is the value in the dacq file. Numbers will be converted
%   to numbers, ambiguous text will remain as strings. Also adds a pixel ratio vector which should match position data
%   hdata = getDACQDATAHEADERS(snames,ftype)
%
%%%%%%%% Inputs
%   snames = the filenames in a cell array
%   ftype = (default = '.pos') the type of file to open, given as an extension e.g. '.pos'
%
%%%%%%%% Outputs
%   hdata = a structure, each row corresponds to an sname, fields are treated as columns
%   e.g. for .pos files the fields include but are not limited to:
%         hdata.trial_time = string corresponding to trial onset
%         hdata.experimenter = text from dacqUSB 'experimenter' text box
%         hdata.comments = text from dacqUSB 'comments' text box
%         hdata.duration = duration of recording in s
%         hdata.num_colours = number of colours tracked
%         hdata.min_x = minimum x value of camera window, max_x, min_y, max_y also exist
%         hdata.window_min_x = minimum x value of tracking window, window_max_x, window_min_y, window_max_y also exist
%         hdata.timebase = frequency of position sampling
%         hdata.pixels_per_metre = pixels per metre value taken from dacqUSB 'pixel ratio' text box
%         hdata.num_pos_samples = the number of position data samples recorded, should be duration*timebase
%
%%%%%%%% Comments
%   21/07/17 created to replace completely stupid headerCheck and key_value functions, to allow variable pixel ratio etc
%   21/07/17 finished .pos and .set conditions
%   24/07/17 added LED colours and bearings from .set file
%   24/07/17 added functionality for all file types, realised that they can all be contained in the same loop
%   24/07/17 regexp has some buggy effects, swicthed to strfind, streamlined conversion of numbers to digits
%
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('ftype','var') || isempty(ftype) || all(isnan(ftype))
    ftype = '.pos';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read files
hdata = struct;
for ss = 1:length(snames) % for every file name
    % open file
    fid = fopen([snames{ss} ftype],'r','ieee-be');
    if fid < 0
        return
    end    

    % find lines with session data
    fseek(fid,0,-1);
    while ~feof(fid) % while we are not at the end of the file
        % get next line of text
        txt = fgetl(fid);
        
        % if this is the file end, skip it
        if strfind(txt,'data_start')
            break
        end
        
        % Otherwise, find where the first space is (occurs right after the variable name)
        sindx = regexp(txt,' ');
        varnow = txt(1:sindx(1)-1);
        valnow = txt(length(varnow)+2:end);
        hdata(ss).(varnow) = valnow;
        
        % deal with variables that are decimals or have extra spaces in them
        dindx = isstrprop(valnow,'digit');
        dindx(strfind(valnow,' ')) = 1;
        if numel(strfind(valnow,'.')) > 1 % if there are multiple decimal places, skip it (can't be a number)
            continue
        else
            dindx(strfind(valnow,'.')) = 1;
        end
        
        % if all the characters are spaces, dots and numbers
        if all(dindx)
            valnow = str2double(valnow); % convert to a number
        end
        hdata(ss).(varnow) = valnow;
    end   
    fclose(fid);  

    % convert text 'timebase' to a number value
    for pp = 1:length(hdata)
        if isfield(hdata(pp),'timebase')
            tbase = hdata(pp).timebase;
            dindx = isstrprop(tbase,'digit');
            tbase = str2double(tbase(dindx));
            hdata(pp).sample_rate_num = tbase;
        end
    end        

    % also add a vector for pixel ratio, if possible, for later functions
    for pp = 1:length(hdata)
        if isfield(hdata(pp),'num_pos_samples') && isfield(hdata(pp),'pixels_per_metre') 
            pnow = ones(hdata(pp).num_pos_samples,1) .* hdata(pp).pixels_per_metre;
            hdata(pp).pixels_per_metre_vec = uint16(pnow);
        end
    end       
end


































