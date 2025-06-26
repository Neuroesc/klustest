function [hdata,fname] = get_set_for_klustest(dataformat,config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getTRODES  find the names of the sessions contributing to this cutfile and the active tetrodes
% Given the combined name of a session this function checks which tetrode has a cut file, opens it and extracts the name of the 
% files (sessions) that contributed towards it. In this process we also get a list of which tetrodes have a .cut file (i.e. were
% cluster cut).
%
% USAGE:
%         [tets,mtets,snames,cutname] = getTRODES(cname,tin)
%
% INPUT:
%         cname - combined name of the session, without an extension, default is 'kwiktint' because thats the default output for kwiktint
%         tin - (optional) tetrodes to check for, default is 1:16
%
% OUTPUT:
%    tets - list of tetrodes with cut files
%    snames - cell array of session names
%    cutname - the name of the cutfile that was used
%
% EXAMPLES:
%
% See also: klustest getcut setdiff

% HISTORY:
% version 1.0.0, Release 24/08/16 Initial release
% version 1.0.1, Release 04/04/18 Commenting and formatting for klustest update
% version 2.0.0, Release 04/04/18 Changed to also search for active tetrodes (combined getTRODES and getCNAMES)
% version 2.0.1, Release 04/04/18 Removed list of inactive tetrodes, this can be generated using setdiff, reduces vargin to just 1
% version 2.0.2, Release 19/04/19 Updated to remove possible duplicates of tetrodes
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
    hdata = table;
    hdata.date = datetime("now",'format','yyyyMMddHHmmss'); % default is now

    switch dataformat
        case {'kwiktint'}
            % d = dir([config.snames{1} '\*.set']);
            % fnames_set = {d.name};
            fname = [config.snames{1} '\' config.snames{1} '.set'];
            h = get_dacq_headers(fname);

            hdata.date = h.trial_date;

        case {'kwikcut','neuralynx'}
            d = dir([config.snames{1} '\*.ntt']);
            fnames_ntt = {d.name};
            n_ntt = length(fnames_ntt);
            d = dir([config.snames{1} '\*.nst']);
            fnames_nst = {d.name};  
            n_nst = length(fnames_nst);        
            fnames = [fnames_ntt; fnames_nst];

            fname = [config.snames{1} '\' fnames{1}];
            h = Nlx2MatSpike(fname,[0 0 0 0 0],1,1,1); 
        
            idx = contains(h,{'TimeCreated'}); % find the time created entry
            if any(idx)
                txt = h{idx}(13:24); % extract the number part of the text
                sindx = regexp(txt,'/');
                txt(sindx) = [];
                hdata.date = txt;     
            end
      
        case {'phy'}
            hdata.date = datestr(now,30);
            fname = 'phy';

    end
    
    
    
    
    
    
    
    
    
    
    







