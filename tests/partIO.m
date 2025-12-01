function part_config = partIO(pname,part_config,oride)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%partIO  loads or saves a part_config for klustest
% Function to load in or save a part_config table (in struct format) in a .mat file
% The reason this is saved in a structure is so that the current part_config is always
% saved in the field named 'part_config' while older versions are saved in different fields
% with fieldnames reflecting the date/time it was overwritten
%
% USAGE:
%         part_config = partIO(pname,part_config,oride)
%
% INPUT:
%         pname - filename of the part_config mat file
%         part_config - a current part_config
%         oride - override setting, if 1 the function will save over an existing part_config file
%
% OUTPUT:
%    part_config - loaded part_config, or the same as input if one was given
%
% EXAMPLES:
%
% See also: klustest save struct

% HISTORY:
% version 1.0.0, Release 16/04/19 Initial release, code to contain loading/saving part_config
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2019 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
    % deal with input variables
    inps = {'pname','part_config','oride'};
    vals = {'klustest\kwiktint\kwiktint_part_config.mat','0','0'};
    reqd = [0 0 0];
    for ff = 1:length(inps)
        if ~exist(inps{ff},'var')
            if reqd(ff)
                error('ERROR: vargin %s missing... exiting',inps{ff});            
            end        
            eval([inps{ff} '=' vals{ff} ';']);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
    if ~exist(pname,'file') % if there is no pre-existing part_config we need to write a new one regardless of the override setting
        if ~istable(part_config)
            error('No part_config given, but no file (%s) exists either, exiting...',pname)
        end
        disp(sprintf('\t...writing new part_config: %s',pname));
        part_data = struct;
        part_data.part_config = part_config;
        save(pname,'part_data'); % save the part_config as a .mat file            

    elseif ~oride % if we do not want to override an existing part_config
        disp(sprintf('\t...loading part_config: %s',pname));
        load(pname,'part_data'); % load the part_config from a .mat file      
        part_config = part_data.part_config;
        
    elseif exist(pname,'file') && oride==1 % if there is a pre-existing part_config but we want to override it
        if ~istable(part_config)
            error('No part_config given, but set to override file (%s), exiting...',pname)
        end        
        disp(sprintf('\t...overwriting part_config: %s',pname));
        load(pname,'part_data'); % load the part_config from a .mat file         
        part_data.(['part_config_' datestr(now,30)]) = part_data.part_config;
        part_data.part_config = part_config;
        save(pname,'part_data'); % save the part_config as a .mat file   
        
    elseif oride==2 % special case, we want to use the current part_config without saving (mainly for debugging)
        if ~istable(part_config)
            error('No part_config given, but set to continue with current part_config, exiting...')
        end           
        disp(sprintf('\t...using current part_config without saving'));

    end









































