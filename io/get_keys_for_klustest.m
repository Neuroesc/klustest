function [keys] = get_keys_for_klustest(dataformat,data_dirs,pos,tstart)
% get_spk_for_klustest load spike waveforms and spike times
% Data loading function for klustest, given a data format type, and some
% settings, loads waveforms and spike times
%
% USAGE
%
% [spk,wav,spk_srate,wavtime] = get_spk_for_klustest(dataformat,data_dirs,tstart)
%
% INPUT
%
% 'dataformat' - String: 'kwiktint'
%
% 'data_dirs' - Data directories, output from get_tets_for_klustest
%
% 'pos' - Table, 'pox','poy','pot','pov','poh','pod','poa'
%       (x, y, time, speed, head direction, displacement direction, angular head
%       velocity)
%
% 'tstart' - Start time of the recordings, given by get_pos_for_klustest
%
% OUTPUT
%
% NOTES
% 1. 
%
% 2. 
% 
% SEE ALSO kwiktint klustest get_dacq_headers

% HISTORY
%
% version 1.0.0, Release 10/04/26 Initial release
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2026 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
%%%%%%%%%%%%%%%% KEY FILES
% run through the data and convert .inp to .key files
% run through the data and actually make the files
keyboard
    for ii = 1:length(snames)
        if revec(ii) % if we need or want to (re)create the key file
            % sort out file names
            [~,nme,~] = fileparts(snames{ii});
            inpname = [pwd '\' nme '.inp'];
            keyname = [pwd '\' nme '.klustest_keys']; 
    
            % read .inp file            
            if ~exist(inpname,'file')
                error(sprintf('ERROR: .inp file %s not found... unable to load',inpname))
            end
            [count,timestamps,type,value] = read_key(inpname);
    
            % cut the data to the manual keypresses only
            mandx = uint32(cell2mat(type(:)))==75;
            if ~any(mandx)
                xin = table;
                xin.timestamps = zeros(0);
                xin.text_value = zeros(0);    
            else      
                value = value(mandx);
                timestamps = timestamps(mandx);
    
                % process the value data to fix errors
                value = cellfun(@transpose,value,'UniformOutput',0); 
                value = cellfun(@num2str,value,'UniformOutput',0);
                value = cellfun(@strtrim,value,'UniformOutput',0);
                value = regexprep(value,{'\W',' '},'');
                value(cellfun(@isempty,value)) = {'#'}; % some keys are not recognised by uint32, replace these with NaNs, usually this happens when we hit a weird key by accident like ` or ¬
                value = value';
    
                xin = table;
                xin.timestamps = cellstr(num2str(timestamps'));
                xin.text_value = cell2mat(value(:));
            end
            % write key file
            writetable(xin,keyname,'FileType','text','Delimiter','tab');
        end    
    end
    












%%%%%%%%%%%%%%%% FIX KEYS
% plot the data and allow the user to fix errors


































