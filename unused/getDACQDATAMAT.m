function tintmat = getDACQDATAMAT(mname,cname,snames,tetrodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	This is a wrapper function for loading of Axona data sets. It aims to create a giant 
%   structure file full of all the session data
%   mtint = getDACQDATA(cname,snames,tetrodes);
%
%%%%%%%% Inputs
%   filepath = the path to the files (usually the working directory)
%   filename = the file names that are to be analysed, if multiple sessions were merged then their names must be given as a cell array
%   cname = this is the name of the output of merging the sessions, I use the two names concatenated in the order they are analysed
%   tetrodes = the tetrodes requested
%
%%%%%%%% Outputs
%   mtint = the structure array of information (details below)
%
%%%%%%%% Comments
%   05/08/16 created from readAllDACQdata
%   06/08/16 added exporting of waveforms for plotting
%   06/08/16 fixed issues with concatenated sessions
%   08/08/16 added use of concatenated filename
%   19/08/16 fixed problem with name, added _multi
%   25/08/16 added tetrodes input, otherwise these functions want to run on all tetrode files, which can be a problem if one is corrupted or not wanted etc
%   26/08/16 changed to getDACQDATA to resolve any ambiguity
%   26/08/16 function now only handles one input file - I will never want to combine files that were not clustered together, I don't see why anyone else would/should either
%
%%%%%%%% mtint details (to the best of my knowledge)
%	mtint.flnmroot = the filename (or date)
%	mtint.filepath = the current directory
%	mtint.header = session settings (.set file?)
%   mtint.channels = individual channel data
%       mtint.channels.dead_channels = which channels are dead (eeg or grounded)
%       mtint.channels.lfp_channels = which channels are used for eeg
%       mtint.channels.grounded_channels = which channels are grounded
%       mtint.channels.lfp_channel_info = a matrix, each row corresponds to an eeg channel, columns: [slot #, channel #, save], see readSET for more info
%       mtint.channels.rec_channel_info = a matrix, each row corresponds to a recording channel, columns: [channel #, mode, filter, gain], see readSET for more info
%	mtint.pos = pos data sub structure
%		mtint.pos.led_pos = the raw positions of the LEDs?
%		mtint.pos.ts = position time
%		mtint.pos.led_pix = 
%		mtint.pos.trial_duration = the duration of the session in seconds
%		mtint.pos.header = more session info, including camera and window settings
%		mtint.pos.xy_pixels = the position of the animal in pixel coordinates?
%		mtint.pos.xy_cm = the position of the animal in coordinates translated to cm?
%		mtint.pos.dir = head direction?
%		mtint.pos.speed = the velocity
%		mtint.pos.jumpyPercen = 
%	mtint.tetrode = spike data sub structure
%		mtint.tetrode.id = electrode #
%		mtint.tetrode.header = session info for each tetrode in turn
%		mtint.tetrode.ts = 
%		mtint.tetrode.ch1 = an array with spikes sorted 1 per row, 50 samples in columns
%		mtint.tetrode.ch2 = an array with spikes sorted 1 per row, 50 samples in columns
%		mtint.tetrode.ch3 = an array with spikes sorted 1 per row, 50 samples in columns
%		mtint.tetrode.ch4 = an array with spikes sorted 1 per row, 50 samples in columns
%		mtint.tetrode.pos_sample = 
%		mtint.tetrode.cut
%	mtint.lfp = local field potential (eeg) sub structure
%		mtint.lfp.lfp = the EEG signal
%		mtint.lfp.Fs = the EEG sampling rate
%		mtint.lfp.header = session information for the eeg data
%   mtint.inps = digital input sub structure (if a .inp file exists)
%       mtint.inps.name = the name of the .inp file
%       mtint.inps.ts = ts;
%       mtint.inps.type = type;
%       mtint.inps.key = key;
%       mtint.inps.header = header;
%
tic1 = tic;
[snames,cname,nsess] = getSNAMES;
mname = [pwd '\3D matfiles\tintmat.mat'];
wname = ['C:\Users\Roddy\Desktop\test\wavmat.mat'];

save_wav = 0; % set to 1 to save a wavmat file, containing all the spike data (probably not necessary as it is quite fast to just read the spike data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect data
disp(sprintf('Running getDACQDATAMAT...'))

% initialise tintmat
tintmat = struct;
wavmat = struct;

% initialise variables
filepath = [pwd '\'];
flnmroot = cname; 
tetrodes = 1:16;

%% Add some metadata to mtint
disp(sprintf('\t...loading metadata'))
sname = [filepath,snames{1},'.set'];

% get the some state info about the recording/eeg channels
[dedch,lfpch,grdch,eblock,cblock] = readSET(sname); % open first .set file to find dead channels - these shouldn't change and Tint can't handle it if they do anyway
tintmat.chn_dead_channels = uint8(dedch);
tintmat.chn_lfp_channels = uint8(lfpch);
tintmat.chn_grounded_channels = uint8(grdch);
tintmat.chn_lfp_channel_info = uint8(eblock);
tintmat.chn_rec_channel_info = uint8(cblock);

% check all required files are present
[fileStruct,tetsAvailable] = listDACQFILES(filepath,snames);

% process headers
tintmat.header = getDACQHeader([filepath,fileStruct(1).flnmroot,'.set'],'set'); % for now use only the first trials .set header

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Load .inp files
disp(sprintf('\t...loading digital inputs'))

% check all inp files are present
inpPresent = zeros(1,length(snames));
for i = 1:length(snames)
    inpPresent(i) = exist([pwd '\' snames{i} '.inp'],'file');
end % for i = 1:length(snames)

% initialise variables
all_in_count = 0;
all_in_tstamps = [];
all_in_type = [];
all_in_value = [];

% read .inp files and accumulate their data
if all(inpPresent)
    for sn = 1:length(snames)
        [count,timestamps,type,value] = read_key([snames{sn} '.inp']);
        all_in_count = all_in_count + count;
        all_in_tstamps = [all_in_tstamps; timestamps(:)];
        all_in_type = [all_in_type; type(:)];
        all_in_value = [all_in_value; value(:)];
    end % for ifile = 1:numel(fileStruct)
end % if all(eegPresent)

% collect data
tintmat.inp_count = single(all_in_count);
tintmat.inp_tstamps = single(all_in_tstamps);
%tintmat.inp_type = all_in_type;
tintmat.inp_value = all_in_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load .pos files
fprintf('\t...loading positions:    ')
idx1 = find(strcmpi('colactive_1',tintmat.header(:,1)));
idx2 = find(strcmpi('colactive_2',tintmat.header(:,1)));
led1 = str2double(char(tintmat.header(idx1,2)));
led2 = str2double(char(tintmat.header(idx2,2)));

% figure out how many LEDs were used
if led1 && led2
    n_leds = 2;
elseif led1
    n_leds = 1;
else
    error('\tERROR: no LEDS tracked... exiting')
end % if led1 & led2

% initialise variables
trial_duration = NaN(1,numel(fileStruct));
total_duration = 0;
trial_samps = NaN(1,numel(fileStruct));
all_led_pos = [];
all_led_pix = [];
all_ts = [];

% retrieve and concatenate contents of .pos files
loopout = loopCOUNT(numel(fileStruct));
for pp = 1:numel(fileStruct)
    current_duration = key_value('duration',getDACQHeader([filepath,fileStruct(pp).flnmroot,'.pos'],'pos'),'num');
    trial_duration(pp) = current_duration;
    
    [led_pos,post,led_pix] = rawpos([filepath,fileStruct(pp).flnmroot,'.pos'],n_leds); % raw led positions
    trial_samps(pp) = length(led_pos(:,1,1));
    
    all_led_pos = [all_led_pos; led_pos];
    all_led_pix = [all_led_pix; led_pix];
    all_ts = [all_ts; post + total_duration];    
    total_duration = total_duration + current_duration;
    loopCOUNT(numel(fileStruct),pp,loopout);
end % for ifile = 1:numel(fileStruct)

% collect data
tintmat.pos_led_pos = int16(all_led_pos);
tintmat.pos_led_pix = uint16(all_led_pix);
tintmat.pos_ts = single(all_ts);
tintmat.pos_total_duration = single(total_duration);
tintmat.pos_trial_duration = single(trial_duration);
tintmat.pos_trial_samps = single(trial_samps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load .cut files
fprintf('\t...loading cluster cuts:')

% load the cut file into the relevant tetrode part of the structure
loopout = loopCOUNT(numel(tetsAvailable));
for tt = 1:numel(tetsAvailable)
    current_tet = tetsAvailable(tt);       
    cutname = [cname '_',num2str(current_tet),'.cut'];

    % if a cutfile exists, then enter its contents into to mtint struct as well - this will allow for analysis of clusters afterwards
    if exist(cutname,'file')
        [clust,exact_text] = getcut([filepath,cutname]); % this reads the .cut file line by line so that it can be entered into the mtint structure  
        numSpikesInCutFile = numel(clust);        
        idx = strfind(exact_text,': '); % should contain two values
        flist = textscan(exact_text(idx(1)+2:idx(2)-8),'%s','delimiter',',');
        fileStruct = sortstruct(fileStruct,'flnmroot',flist{1});
    else
        clust = [];
    end % if exist([filepath,fileStruct(1).flnmroot,'_',num2str(current_tet),'.cut'],'file') 
    
    % collect data
    tintmat.(['tet_' num2str(current_tet) '_cut']) = uint8(clust);
    tintmat.(['tet_' num2str(current_tet) '_nspikes']) = numSpikesInCutFile;    
    loopCOUNT(numel(tetsAvailable),tt,loopout);
end % for itet = 1:numel(mtint.tetrode)

for i = 1:numel(fileStruct)
    filelist(i) = cellstr(fileStruct(i).flnmroot);
end % for i = 1:numel(fileStruct)

% use the first file as a header
tintmat.flnmroot = filelist;
tintmat.filepath = filepath;
tintmat.pos_header = headerCheck(filepath,filelist,'pos');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Load tetrodes
fprintf('\t...loading spikes:        ')
durations = zeros(1,numel(fileStruct)+1);
loopout = loopCOUNT(numel(tetsAvailable));
for tt = 1:numel(tetsAvailable)
    current_tet = tetsAvailable(tt);     
    nspikes_cut = tintmat.(['tet_' num2str(current_tet) '_nspikes']);
    
    % load feature file (this is faster than reading the spike files
    fetfile = [cname '.fet.' num2str(current_tet)];
    fid = fopen(fetfile,'r');
    fets = sscanf(fgetl(fid),'%d');
    nfets = (fets-1)/4;
    fdat = fscanf(fid,'%f',[fets Inf]);
    fdata = fdat';            
    tet_ts = fdata(:,end);

    % work out the corresponding position samples
    tet_pos_sample = ceil(tet_ts * 50);

    % collect data
    tintmat.(['tet_' num2str(current_tet) '_ts']) = single(tet_ts);
    tintmat.(['tet_' num2str(current_tet) '_pos_sample']) = single(tet_pos_sample);
    if save_wav % if we want to save waveform information
        wavmat.(['tet_' num2str(current_tet) '_ts']) = single(tet_ts);
        wavmat.(['tet_' num2str(tt) '_ch1']) = tet_ch1;
        wavmat.(['tet_' num2str(tt) '_ch2']) = tet_ch2;
        wavmat.(['tet_' num2str(tt) '_ch3']) = tet_ch3;
        wavmat.(['tet_' num2str(tt) '_ch4']) = tet_ch4;
    end % if save_wav % if we want to save waveform information

    % check if getspikes and the cut file agree
    nspikes_tet = numel(tintmat.(['tet_' num2str(current_tet) '_ts']));
    if ~isequal(nspikes_tet,nspikes_cut)
        error('ERROR: nspikes in tetrode (%d) and nspikes in cut file (%d) do not match foir tetrode %d... exiting',nspikes_tet,nspikes_cut,num2str(current_tet));
    end % if ~isequal(nspikes_tet,nspikes_clust)
    loopCOUNT(numel(tetsAvailable),tt,loopout);
end % for itet = 1:numel(tetsAvailable)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Load .eeg files
fprintf('\t...loading lfp:            ')

% list all eeg files
eegfilelist = dir([filepath,snames{1},'.eeg*']);

% read them in and assign to structure
loopout = loopCOUNT(numel(eegfilelist));
for ifile = 1:numel(eegfilelist)
    enow = eegfilelist(ifile).name;
    [~,~,ext] = fileparts(enow);

    % check all eeg files are present
    eegPresent = zeros(1,length(snames));
    for i = 1:length(snames)
        eegPresent(i) = exist([pwd '\' snames{i} ext],'file');
    end % for i = 1:numel(fileStruct)

    if all(eegPresent)
        for sn = 1:length(snames)
            [EEG,Fs] = geteeg([filepath,snames{sn},ext]);
            tintmat.(['lfp_' num2str(ifile) '_eeg']) = ext; % add the extension of the eeg file    
            tintmat.(['lfp_' num2str(ifile) '_lfp']) = int8(EEG); % concatenate eeg
            tintmat.(['lfp_' num2str(ifile) '_Fs']) = Fs; % add the sampling rate
        end % for ifile = 1:numel(fileStruct)
    end % if all(eegPresent)
    loopCOUNT(numel(eegfilelist),ifile,loopout);
end % for ifile = 1:numel(eegfilelist)





%     % Initialise variables
%     tet_ts = int8(zeros(nspikes_cut,1)); 
%     if save_wav % if we want to save waveform information
%         tet_ch1 = int8(zeros(nspikes_cut,50));
%         tet_ch2 = int8(zeros(nspikes_cut,50));
%         tet_ch3 = int8(zeros(nspikes_cut,50));
%         tet_ch4 = int8(zeros(nspikes_cut,50)); 
%     end % if save_wav

    % retrieve and concatenate contents of .tet files
%     read_spikes = 1;
%     for ifile = 1:numel(fileStruct)
%         durations(ifile+1) = key_value('duration',getDACQHeader([filepath,fileStruct(ifile).flnmroot,'.',num2str(current_tet)],'tet'),'num');
%         duration = sum(durations(1:ifile));
%         if save_wav % if we want to save waveform information
%             [ts,ch1,ch2,ch3,ch4] = getspikes([filepath,fileStruct(ifile).flnmroot,'.',num2str(current_tet)]);            
%             ts = ts + duration;
%             tet_ts(read_spikes:read_spikes+numel(ts)-1,:) = ts;
%             tet_ch1(read_spikes:read_spikes+numel(ts)-1,:) = ch1;
%             tet_ch2(read_spikes:read_spikes+numel(ts)-1,:) = ch2;
%             tet_ch3(read_spikes:read_spikes+numel(ts)-1,:) = ch3;
%             tet_ch4(read_spikes:read_spikes+numel(ts)-1,:) = ch4;        
%             read_spikes = read_spikes+numel(ts);
%         else
%             fetfile = [cname '.fet.' num2str(current_tet)];
%             fid = fopen(fetfile,'r');
%             fets = sscanf(fgetl(fid),'%d');
%             nfets = (fets-1)/4;
%             fdat = fscanf(fid,'%f',[fets Inf]);
%             fdata = fdat';            
%             ts = fdata(:,end);
%             
%             
%             
%             [ts,~,~,~,~] = getspikes([filepath,fileStruct(ifile).flnmroot,'.',num2str(current_tet)],1);            
%             ts = ts + duration;
%             tet_ts(read_spikes:read_spikes+numel(ts)-1,:) = ts;            
%         end % if save_wav % if we want to save waveform information
%     end % for ifile = 1:numel(fileStruct)

