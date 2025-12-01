function part_config = partLOAD(pname,part_config,oride)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION  short descr.
% long descr.
%
% USAGE:
%         [out] = template(in,in2)
%
% INPUT:
%         in - input 1
%         in2 - input 2
%
% OUTPUT:
%    p - output
%
% EXAMPLES:
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2019 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with input variables
inps = {'pname','oride'};
vals = {'klustest\kwiktint\kwiktint_part_config.mat','0'};
reqd = [0 0];
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
    if ~exist(pconfig_name,'file')
        oride = 1;
    end

    % load or save part_config
    switch oride
        case {0}
            disp(sprintf('\t...loading part_config'));
            load(pconfig_name,'part_config'); % load the part_config from a .mat file            

        case {1}
            disp(sprintf('\t...writing new part_config'));
            save(pconfig_name,'part_config'); % save the part_config as a .mat file            

        case {2}
            disp(sprintf('\t...using current part_config without saving'));
    end





















































