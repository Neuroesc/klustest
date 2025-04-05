function [d,e,g,eblock,cblock,lfpnames] = readSET(setname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function reads a .set file, mainly to ascertain which channels were used for EEG or grounded
%   The reason is that grounded channels passed to Tint for klustakwiking will cause it to find only 
%   one cluster every time. Strangely, EEG channels are recorded as being grounded, so the function also
%   outputs the uniquely dead channels.
%   [d,e,g] = readSET(setname)
%
%%%%%%%% Inputs
%   setname = the full name of the setfile (can be a directory address if necessary)
%
%%%%%%%% Outputs
%   dedch = a vector containing all the dead channels (i.e. grounded or EEG) with no repetition
%   eegch = a vector containing the channels used for EEG
%   grdch = a vector containing the channels that were grounded
%
%   eblock = a matrix, each row corresponds to an eeg channel, columns: [slot #, channel #, save]
%           slots range from 1:64 and are the rows in the dacqUSB menu for eeg channels
%           channels correspond to recording channels 1:64
%           save is 1 if the eeg has been ticked (to save) in dacqUSB, 0 otherwise
%   cblock = a matrix, each row corresponds to a recording channel, columns: [channel #, mode, filter, gain]
%           channels correspond to recording channels 1:64
%           mode: 0 = signal, 0 = ref, 2 = -signal, 3 = -ref, 4 = sig-ref, 5 = ref-sig, 6 = grounded
%           filter: 0 = direct, 2 = direct+notch, 2 = highpass, 3 = lowpass, 4 = lowpass+notch
%           gain is the gain of the channel
%   lfpnames = a cell array, each row represents an eeg channel, columns: [eeg name, egf name, actual channel recorded]
%
%%%%%%%% Comments
%   18/08/16 created to handle channel inputs to Tint command line interface
%   24/08/16 fixed a bug
%   25/08/16 majorly revised to handle the weird dacqUSB menu system, added eblock and cblock output
%   25/08/16 added checks for common mistakes in .set files
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('setname','var') || isempty(setname)
    error('ERROR: readSET is missing .set name... exiting')
end % if ~exist('setname','var') || isempty(setname)

%% Open file
fID = fopen(setname,'r');
[~,fname,~] = fileparts(setname);

sdata = textscan(fID,'%s','Delimiter','\r');

eblock = NaN(128,3); % will contain eeg data, [slot, channel, save]
cblock = NaN(128,6); % will contain channel data, [channel, mode, filter, gain]
for ll = 1:length(sdata{1,1}) % for every line
    lnow = sdata{1,1}{ll,1}; % get the line
    sindx = strfind(lnow,' '); % index of spaces in text
    uindx = strfind(lnow,'_');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %% deal with eeg channels
    if strfind(lnow,'EEG_ch_'); % if this line descibes eeg
        if isempty(strfind(lnow,'BPF')); % if it is not about the BPF mode
            ch_now = str2double(lnow(uindx(end)+1:sindx(end)));
            if strfind(lnow,'save'); % if it is the 'save' setting
                eblock(ch_now,1) = ch_now;
                eblock(ch_now,3) = str2double(lnow(sindx(end)+1:end));
            else % if it is the 'channel' setting
                eblock(ch_now,1) = ch_now;
                eblock(ch_now,2) = str2double(lnow(sindx(end)+1:end));
            end
        end
    end 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% deal with channel modes
    if strfind(lnow,'mode_ch_'); % if this line descibes channel modes
        if isempty(strfind(lnow,'disp')); % if it is not about the display mode
            ch_now = str2double(lnow(uindx(end)+1:sindx(end)))+1;
            cblock(ch_now,1) = ch_now;
            cblock(ch_now,2) = str2double(lnow(sindx(end)+1:end));
        end 
    end 
    
    %% deal with channel filters
    if strfind(lnow,'filter_ch_'); % if this line descibes channel modes
        ch_now = str2double(lnow(uindx(end)+1:sindx(end)))+1;
        cblock(ch_now,3) = str2double(lnow(sindx(end)+1:end));
    end   
    
    %% deal with channel gain
    if strfind(lnow,'gain_ch_'); % if this line descibes channel modes
        ch_now = str2double(lnow(uindx(end)+1:sindx(end)))+1;
        cblock(ch_now,4) = str2double(lnow(sindx(end)+1:end));
    end    
    
    %% deal with channel signal
    if strfind(lnow,'a_in_ch_'); % if this line descibes what channel this slot is actually showing
        ch_now = str2double(lnow(uindx(end)+1:sindx(end)))+1;
        cblock(ch_now,5) = str2double(lnow(sindx(end)+1:end))+1;
    end       
    
    %% deal with channel reference
    if strfind(lnow,'b_in_ch_'); % if this line descibes what channel this slot is referenced to
        ch_now = str2double(lnow(uindx(end)+1:sindx(end)))+1;
        cblock(ch_now,6) = str2double(lnow(sindx(end)+1:end))+1;
    end  
end 

e = eblock(eblock(:,3)==1,2); % find the eegs set to save (including 0s)
g = cblock(cblock(:,2)==6,1); % find the channels set to grounded
% eeg slots can be set to save eeg from channel 0 or below, in this case it just saves an empty eeg and leaves all the channels alone, so we don't care about that
d = unique([e(:);g(:)]);
d = d(d > 0); % 0 can't be a dead channel so ignore it

%% Work out the 'actual' eegs, i.e. if a channel is set to ref then it records the reference not the slot
eindx = find(eblock(:,3) == 1); % find which eeg slots were set to save
act_eeg = cell(length(e),3);
for ee = 1:length(eindx)
    enow = eindx(ee); % the eeg slot we are looking at
    cnow = eblock(enow,2); % the channel this slot was set to
    if ee == 1
        act_eeg{ee,1} = [fname '.eeg'];
        act_eeg{ee,2} = [fname '.egf'];          
    else
        act_eeg{ee,1} = [fname '.eeg' num2str(ee)];
        act_eeg{ee,2} = [fname '.egf' num2str(ee)];  
    end % if ee == 1

    if cnow <= 0
        act_eeg{ee,3} = 0; % the actual channel was the reference    
    else
        if cblock(cnow,2) == 1 || cblock(cnow,2) == 3 % if this channel's mode was set to ref or -ref      
            act_eeg{ee,3} = cblock(cnow,6); % the actual channel was the reference
        else
            act_eeg{ee,3} = cblock(cnow,5); % the actual channel was just the channel  
        end % if cblock(cnow,2) == 1 || cblock(cnow,2) == 3 % if this channel's mode was set to ref or -ref        
    end % if cnow <= 0 
end % for ee = 1:length(eindx)
lfpnames = act_eeg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do some checks
windx = [find(cblock(:,2) == 5 & cblock(:,3) == 3) find(cblock(:,2) == 5 & cblock(:,3) == 4) find(cblock(:,2) == 4 & cblock(:,3) == 3) find(cblock(:,2) == 4 & cblock(:,3) == 4)]; % find channels set to lowpass or lowpass+notch, but also set to ref-sig or sig-ref
if numel(windx) > 0
    disp(sprintf('\tFIX THIS: channel %s is set to lowpass but is also set to subtract reference',mat2str(windx')))
end % if find(cblock(:,2) == 3 || cblock(:,2) == 4 && cblock(:,3) == 5 || cblock(:,3) == 4) 

windx = find(eblock(:,3) == 1 & eblock(:,2) == 0); % find channels set to save eeg, but where the eeg is zero
if numel(windx) > 0
    disp(sprintf('\tFIX THIS: eeg slot %s is set to save eeg, but the channel is set to 0',mat2str(windx')))
end % if find(cblock(:,2) == 3 || cblock(:,2) == 4 && cblock(:,3) == 5 || cblock(:,3) == 4) 











