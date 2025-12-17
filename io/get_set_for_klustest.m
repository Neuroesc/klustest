function [hdata,fname] = get_set_for_klustest(dataformat,config)
% get_set_for_klustest load session info from .set files or equivalent
% Data loading function for klustest, given a data format type, and some
% settings, finds .set files (Tint format) or .ntt files (Neuralynx format) and
% extracts some headers/information that are useful for analysis,
%
% USAGE
%
% [hdata,fname] = get_set_for_klustest(dataformat,config)
%
% INPUT
%
% 'dataformat' - String: 'kwikcut','neuralynx','kwiktint','phy' (experimental)
%
% 'config' - Structure, configuration file from klustest, contains the
%           input/output filename to check for tetrodes etc
%
% OUTPUT
%
% 'hdata' - Table, session data, headers
%
% 'fname' - Filename of the file used to build hdata
%
% NOTES
% 1. Phy data format is experimental
%
% 2. Function is not really intended to be used without klustest
% 
% SEE ALSO kwiktint klustest Nlx2MatSpike

% HISTORY
%
% version 1.0.0, Release 24/08/16 Initial release modified from getSET
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
    hdata = table;
    hdata.date = datetime("now",'format','yyyyMMddHHmmss'); % default is now

    switch lower(dataformat)
        case {'kwiktint'}
            % d = dir([config.snames{1} '\*.set']);
            % fnames_set = {d.name};
            fname = [config.snames{1} '.set'];
            h = get_dacq_headers(fname);

            hdata.date = h.trial_date;

        case {'kwikcut','neuralynx'}
            d = dir( fullfile(config.snames{1},'*.ntt') );
            fnames_ntt = {d.name};
            n_ntt = length(fnames_ntt);
            d = dir( fullfile(config.snames{1},'*.nst') );
            fnames_nst = {d.name};  
            n_nst = length(fnames_nst);        
            fnames = [fnames_ntt; fnames_nst];

            fname = fullfile(config.snames{1},fnames{1});
            h = Nlx2MatSpike(fname,[0 0 0 0 0],1,1,1); 
        
            idx = contains(h,{'TimeCreated'}); % find the time created entry
            if any(idx)
                txt = h{idx}(13:24); % extract the number part of the text
                sindx = regexp(txt,filesep);
                txt(sindx) = [];
                hdata.date = txt;     
            end
      
        case {'phy'}
            hdata.date = datestr(now,30);
            fname = 'phy';

    end
    
    
    
    
    
    
    
    
    
    
    







