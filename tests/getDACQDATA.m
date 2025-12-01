function mtint = getDACQDATA(config,snames,tetrodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	This is a wrapper function for loading of Axona data sets. It aims to create a giant 
%   structure file full of all the session data
%   mtint = getDACQDATA(cname,snames,tetrodes);
%
%%%%%%%% Inputs
%   cname = this is the name of the output of merging the sessions, i.e. 'kwiktint'
%   snames = the different filenames used to generate the .cut file
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
%   01/09/16 added digital inputs to mtint
%   10/09/16 added LFP to mtint
%   20/09/16 function will now get spike times from the last column of .fet files if they exist (this is faster than reading spike files)
%   21/09/16 added progress counter to spike loading
%   01/11/16 reduced mtint size by converting most things to single, or uint where possible, mtint should now be about 50mB instead of gigabytes
%   10/04/17 added warning for when the .cut file is not found 
%       13/04/16 changed this to an error as it tends to mess everything up when this is the case
%   13/04/17 added an error message for when the .cut file doesnt match the parent filenames
%   07/06/17 fixed error where the function always tries to use all tetrodes (added tetrodes as an input)
%   24/07/17 added compatibility with variable pixel ratio. now using getDACQDATA instead of key_value, also cleaned up headers
%
%%%%%%%% mtint details (additional fields are added during post-processing)
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
%		mtint.pos.led_pix = the size of the light in number of pixels 
%		mtint.pos.trial_duration = the duration of the session in seconds
%		mtint.pos.header = more session info, including camera and window settings
%	mtint.tetrode = spike data sub structure
%		mtint.tetrode.id = electrode #
%		mtint.tetrode.ts = time of every spike on this tetrode
%		mtint.tetrode.cut = cut vector corresponding to spikes on this tetrode
%	mtint.lfp = local field potential (eeg) sub structure
%		mtint.lfp.lfp = the EEG signal
%		mtint.lfp.Fs = the EEG sampling rate
%		mtint.lfp.header = session information for the eeg data
%   mtint.inps = digital input sub structure (if a .inp file exists)
%       mtint.inps.count = the number of timestamps/keystrokes
%       mtint.inps.tstamps = time stamps of key strokes
%       mtint.inps.value = the values of the keystrokes
%
%   © Roddy Grieves: rmgrieves@gmail.com (not sure who wrote the original readalldacqdata though)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
    filepath = [pwd '\'];
    cname = config.cname;
    flnmroot = cname;
    mtint = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(sprintf('\t...running getDACQDATA'))
    disp(sprintf('\t...getting metadata'))

%% Add some metadata to mtint
    mtint.flnmroot = flnmroot;
    mtint.filepath = filepath;
    [~,nme,~] = fileparts(snames{1});
    sname = [filepath,nme,'.set'];
    mtint.header = getDACQDATAHEADERS(snames,'.set');

% in some cases I want a numeric date for the recording, so we have to convert the text string given by dacq
    dstring = mtint.header(1).trial_date; % get the text date
    dcells = strsplit(dstring); % split by spaces
    dstring2 = [dcells{1,2} dcells{1,3} dcells{1,4}]; % recombine using the interesting parts
    formatIn = 'ddmmmyyyy'; % the format it is in now
    dnum = datenum(dstring2,formatIn); % convert to numeric time format
    dstring3 = datestr(dnum,'ddmmyyyy'); % convert to a new format with our specification
    mtint.header_date = dstring3; % save this to mtint

% get the some state info about the recording/eeg channels
    [dedch,lfpch,grdch,eblock,cblock] = readSET(sname); % open first .set file to find dead channels - these shouldn't change and Tint can't handle it if they do anyway
    mtint.channels.dead_channels = dedch;
    mtint.channels.lfp_channels = lfpch;
    mtint.channels.grounded_channels = grdch;
    mtint.channels.lfp_channel_info = eblock;
    mtint.channels.rec_channel_info = cblock;

% check all required files are present
    [fileStruct,tetsAvailable] = listDACQFILES(filepath,snames);
    if exist('tetrodes','var') || ~isempty(tetrodes) || ~all(isnan(tetrodes))
        tetsAvailable = tetrodes;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Load .inp files
    disp(sprintf('\t...loading digital inputs'))

% read them in and assign to structure
    inps = struct;
    inps.count = 0; 
    inps.tstamps = []; 
    inps.value = [];
    inps.text_value = [];
    inps.index = [];
    inps.spec = [];

% check all inp files are present
    inpPresent = zeros(1,length(snames));
    for i = 1:length(snames)
        inpPresent(i) = exist([pwd '\' snames{i} '.inp'],'file');
    end

% read .inp files and accumulate their data
    ttime = 0;
    if all(inpPresent)
        for pp = 1:numel(fileStruct)      
            % load raw digital inputs from .inp files
            inpname = [snames{pp} '.inp'];
            if ~exist(inpname,'file')
                error(sprintf('ERROR: .inp file %s not found... unable to load',inpname))
            end            
            [count,timestamps,type,value] = read_key(inpname);

            inps.count = inps.count + length(count);
            timestamps = timestamps + ttime;
            inps.tstamps = [inps.tstamps; single(timestamps(:))];
            inps.value = [inps.value; value(:)];
            inps.index = [inps.index; ones(size(timestamps(:))).*pp];
            timenow = mtint.header(pp).duration;
            ttime = ttime + timenow;
        end 
    end 
    mtint.inps = inps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load .pos files
    disp(sprintf('\t...loading positions'))
    led1 = mtint.header(1).colactive_1;
    led2 = mtint.header(1).colactive_2;

% figure out how many LEDs were used
    if led1 && led2
        n_leds = 2;
    elseif led1
        n_leds = 1;
    else
        error('\tERROR: no LEDS tracked... exiting')
    end

% retrieve and concatenate contents of .pos files
    pos.header = getDACQDATAHEADERS(snames,'.pos');
    pos.led_pos = [];
    pos.led_pix = [];
    pos.ts = [];

    trial_duration = NaN(1,numel(fileStruct));
    total_duration = 0;
    trial_samps = NaN(1,numel(fileStruct));
    loopout = looper(numel(fileStruct));
    for pp = 1:numel(fileStruct)
        trial_duration(pp) = mtint.header(pp).duration;    
%         [led_pos,post,led_pix] = rawpos([filepath,fileStruct(pp).flnmroot,'.pos']); % raw led positions
        [led_pos,post,led_pix] = read_rawpos([filepath,fileStruct(pp).flnmroot,'.pos'],1);
        trial_samps(pp) = length(led_pos(:,1,1));

        pos.led_pos = [pos.led_pos; single(led_pos)];
        pos.led_pix = [pos.led_pix; uint8(led_pix)];
        ts = post + total_duration;
        pos.ts = [pos.ts; single(ts)];
        pos.trial_duration = trial_duration; % in seconds...
        pos.trial_samples = trial_samps; % number of samples per session (important for rectifying stray position data in 3D tracking)

        total_duration = total_duration + trial_duration(pp);
        pos.total_duration = total_duration;
        loopout = looper(loopout);
    end
    mtint.pos = pos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load .cut files
    disp(sprintf('\t...loading cluster cuts'))
    mtint.tetrode = struct('id',{},'ts',{},'cut',{});
    
% load the cut file into the relevant tetrode part of the structure
    loopout = looper(length(tetsAvailable));
    for tt = 1:length(tetsAvailable) % for every possible tetrode  
        current_tet = tetsAvailable(tt);
        cutname = [cname '_',num2str(current_tet),'.cut'];

        % if a cutfile exists, then enter its contents into to mtint struct as well - this will allow for analysis of clusters afterwards
        if exist(cutname,'file')
            [clust,exact_text] = getcut([filepath,cutname]); % this reads the .cut file line by line so that it can be entered into the mtint structure  
            numSpikesInCutFile = numel(clust);
            idx = strfind(exact_text,': '); % should contain two values
            flist = textscan(exact_text(idx(1)+2:idx(2)-8),'%s','delimiter',',');

            for ff = 1:length(flist{1})
                if ~any(strcmp(flist{1}{ff},{fileStruct.flnmroot}.'))
                    error('\tERROR: filename (%s) does not seem to match any filenames in .cut file (%s)...',flist{1}{ff},strjoin({fileStruct.flnmroot}.',', '))
                end
            end

            fileStruct = sortstruct(fileStruct,'flnmroot',flist{1});
        else
            error('ERROR: no cut file found for tetrode %d... check filename (%s) is correct or define input tetrodes',current_tet,cutname);
        end
        mtint.tetrode(current_tet).id = current_tet;
        mtint.tetrode(current_tet).cut = uint8(clust);
        mtint.tetrode(current_tet).nspike_cut = numSpikesInCutFile; 
        loopout = looper(loopout);
    end

    for i = 1:numel(fileStruct)
        filelist(i) = cellstr(fileStruct(i).flnmroot);
    end

% use the first file as a header
    mtint.flnmroot = filelist;
    mtint.filepath = filepath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Load tetrodes
    disp(sprintf('\t...loading spikes'));
    durations = zeros(1,numel(fileStruct)+1);
    msg = 0;
    loopout = looper(length(tetsAvailable));
    for tt = 1:length(tetsAvailable) % for every possible tetrode  
        current_tet = tetsAvailable(tt);  
        nspikes_cut = mtint.tetrode(current_tet).nspike_cut;

        if isnan(nspikes_cut)
            mtint.tetrode(current_tet).ts = NaN;
            continue
        end

        fetfile = ([filepath,cname,'.fet.',num2str(current_tet)]);
        if exist(fetfile,'file') % if we can use a .fet file use this because it is generally faster (time of each spike is last column)
            msg = 1;
            fid = fopen(fetfile,'r');
            fets = sscanf(fgetl(fid),'%d');
            nfets = (fets-1)/4;
            fdat = fscanf(fid,'%f',[fets Inf]);
            fdata = fdat';            
            tet_ts = fdata(:,end); 
            mtint.tetrode(current_tet).ts = single(tet_ts);       
        else
            all_ts = single(zeros(nspikes_cut,1));
            scount = 1;
            for ifile = 1:numel(fileStruct)
                tetfile = ([filepath,fileStruct(ifile).flnmroot,'.',num2str(current_tet)]);

                durations(ifile+1) = mtint.header(ifile).duration;    
                duration = sum(durations(1:ifile));
                mtint.tetrode(current_tet).id = current_tet;
                [ts,~,~,~,~] = getspikes([filepath,fileStruct(ifile).flnmroot,'.',num2str(current_tet)]);            
                ts = ts + duration;
                all_ts(scount:scount+numel(ts)-1) = single(ts);
                scount = scount+numel(ts);
            end 
            mtint.tetrode(current_tet).ts = all_ts;
        end 

        nspikes_tet = numel(mtint.tetrode(current_tet).ts);
        if ~isequal(nspikes_tet,nspikes_cut)
            error('ERROR: tetrode %d number of spikes in cut file (%d) is different to number of spikes in data file (%d)... exiting',current_tet,nspikes_cut,nspikes_tet);
        end 
        loopout = looper(loopout);
    end 
    if msg
        disp(sprintf('\b (used .fet files)'))
    end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Load .eeg files
    disp(sprintf('\t...loading lfp'))
    eegfilelist = dir([filepath,snames{1},'.eeg*']); % list all eeg files
    lfp = struct('eeg',{},'lfp',{},'Fs',{}); % read them in and assign to structure
    loopout = looper(numel(eegfilelist));
    for ifile = 1:numel(eegfilelist)
        enow = eegfilelist(ifile).name;
        [~,~,ext] = fileparts(enow);
        lfp(ifile).dummy = 1;

        % check all eeg files are present
        eegPresent = zeros(1,length(snames));
        for i = 1:length(snames)
            eegPresent(i) = exist([pwd '\' snames{i} ext],'file');
            if eegPresent(i)
                s = dir([pwd '\' snames{i} ext]);
                if s.bytes < 2^10
                    eegPresent(i) = false;
                end
            end
        end 

        if all(eegPresent)
            for sn = 1:length(snames)
                [lfpv,~,Fs] = getLFPV([filepath,snames{sn},ext]); % load the lfp in microvolts
                lfp(ifile).eeg = ext; % add the extension of the eeg file         
                lfp(ifile).lfp = [lfp(ifile).lfp; single(lfpv)]; % concatenate eeg
                lfp(ifile).Fs = Fs; % add the sampling rate
            end 
        end 
        loopout = looper(loopout);
    end 
    mtint.lfp = lfp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Calculate cluster quality
    disp(sprintf('\t...getting cluster quality'))
    loopout = looper(length(tetsAvailable));
    for tt = 1:length(tetsAvailable) % for every possible tetrode 
        tet = tetsAvailable(tt);  
        fetname = [cname '.fet.' num2str(tet)];
        cutname = [cname '_'  num2str(tet) '.cut'];
        fetdata = clusterQUALITY_v2(fetname,cutname); 
        mtint.clu_quality(tet) = fetdata;
        loopout = looper(loopout);
    end



























