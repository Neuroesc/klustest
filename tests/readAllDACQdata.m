function mtint = readAllDACQdata ( varargin )
%	mtint = readAllDACQdata (filepath,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	% Wrapper function for loading of Axona data sets. 
%	% Creates a giant struct file full of all the session data
%
%	mtint.flnmroot = the filename (or date)
%	mtint.filepath = the current directory
%	mtint.header = session settings (.set file?)
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
%	mtint.EEG = EEG sub structure
%		mtint.EEG.EEG = the EEG signal
%		mtint.EEG.Fs = the EEG sampling rate
%	mtint.eeg = another EEG sub structure
%		mtint.eeg.header = session information for the eeg data
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
	[filename,filepath] = uigetfile('*.set','Select the .set file...',...
        'MultiSelect','on');
elseif nargin == 2
    	filepath = varargin{1};
    	filename = varargin{2};
end % if nargin == 0
if ~strcmp(filepath(end),filesep)
    	filepath(end+1) = filesep;
end % if ~strcmp(filepath(end),filesep)
% check all required files are present
[fileStruct,tetsAvailable] = list_and_check_DACQFiles( filepath,cellstr(filename) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(fileStruct) > 1
% process headers
% headerCheck() - concatenates the values contained in the header files
% of multiple files e.g. sums durations, takes min x values for tracker
% min x etc etc
    	posHeaders = headerCheck(filepath,cellstr(fileStruct(1).flnmroot),'pos');
    	pos = struct('led_pos',{},'ts',{},'led_pix',{},'header',{});
    	pos(1).header = posHeaders;
% for now use only the first trials .set header
    	mtint.header = getDACQHeader ( [filepath,fileStruct(1).flnmroot,'.set'], 'set' );
% load pos
    	idx1 = find(strcmpi('colactive_1',mtint.header(:,1)));
    	idx2 = find(strcmpi('colactive_2',mtint.header(:,1)));
    	led1 = str2double(char(mtint.header(idx1,2)));
    	led2 = str2double(char(mtint.header(idx2,2)));
    	if led1 & led2
       		n_leds = 2;
    	elseif led1
        	n_leds = 1;
    	end % if led1 & led2

    	for ifile = 1:numel(fileStruct)
        	if ifile == 1
			duration = 0;
			TrialDuration = key_value('duration',getDACQHeader([filepath,fileStruct(ifile).flnmroot,'.pos'],'pos'),'num');
			trial_duration(ifile)=[TrialDuration];clear TrialDuration
        	else
			current_duration = key_value('duration',getDACQHeader([filepath,fileStruct(ifile).flnmroot,'.pos'],'pos'),'num');
			duration = duration + current_duration;
			trial_duration(ifile)=[current_duration];
        	end % if ifile == 1
		[led_pos,post,led_pix] = rawpos([filepath,fileStruct(ifile).flnmroot,'.pos'],n_leds); % raw led positions
		pos.led_pos = [pos.led_pos;led_pos];
		pos.led_pix = [pos.led_pix;led_pix];
		ts = post + duration;
		pos.ts = [pos.ts;ts];
		pos.trial_duration = trial_duration;%in seconds...
    	end % for ifile = 1:numel(fileStruct)

	clear duration current_duration
	mtint.pos = pos;
	% load tetrodes
	tetrode = struct('id',{},'ts',{},'pos_sample',{},'cut',{},'duration',{});
	durations = zeros(1,numel(fileStruct)+1);

	for itet = 1:numel(tetsAvailable)
		for ifile = 1:numel(fileStruct)
			durations(ifile+1) = key_value('duration',getDACQHeader([filepath,fileStruct(ifile).flnmroot,'.',num2str(tetsAvailable(itet))],'tet'),'num');
			duration = sum(durations(1:ifile));
			tetrode(itet).id = tetsAvailable(itet);
			[ts,~,~,~,~] = getspikes([filepath,fileStruct(ifile).flnmroot,'.',num2str(tetsAvailable(itet))]);
			ts = ts + duration;
			tetrode(itet).ts = [tetrode(itet).ts;ts];
			tetrode(itet).cut = [];
		end % for ifile = 1:numel(fileStruct)
		tetrode(itet).duration = sum(durations);
		clear duration current_duration
		tetrode(itet).pos_sample = ceil(tetrode(itet).ts * 50);
	end % for itet = 1:numel(tetsAvailable)

	mtint.tetrode = tetrode;
	% load eeg
	eeg_out = struct('EEG',{},'Fs',{},'dummy',{});
	eeg_out(1).dummy = 1;
	% check all eeg files are present first
	for i = 1:numel(fileStruct)
		eegPresent = fileStruct(i).eegFiles;
	end % for i = 1:numel(fileStruct)
	if all(eegPresent)
		for ifile = 1:numel(fileStruct)
			[EEG,Fs] = geteeg([filepath,fileStruct(ifile).flnmroot,'.eeg']);
			eeg_out.EEG = [eeg_out.EEG;EEG];
			eeg_out.Fs = [eeg_out.Fs;Fs];
		end % for ifile = 1:numel(fileStruct)
	end % if all(eegPresent)

	eeg_out = rmfield(eeg_out,'dummy');
	mtint.EEG = eeg_out;
	% load the cut file into the relevant tetrode part of the structure
	% for multiple files the cut files in tint classic are named after the
	% first trial loaded in the set. the same convention is used here
	for itet = 1:numel(mtint.tetrode)
		current_tet = mtint.tetrode(itet).id;
		if exist([filepath,fileStruct(1).flnmroot,'_',num2str(current_tet),'.cut'],'file')
			[clust,exact_text] = getcut([filepath,fileStruct(1).flnmroot,'_',num2str(current_tet),'.cut']);
			% check for the correct length of the cut file and that the files
			% specified in the cut file match with the order files were picked
			% here
			numSpikesInCutFile = numel(clust);
			numSpikesInTetrode = numel(mtint.tetrode(itet).ts);
			if numSpikesInCutFile ~= numSpikesInTetrode
				warning('SpikeCountMismatch:cutFile','There aren''t the same number of spikes loaded across trials as there are in the cut file');
				button = questdlg('Do you want to load another cut file?');
				if strcmpi(button,'yes')
					[flnmroot,filepath] = uigetfile('*.cut','Select the .cut file...','MultiSelect','off');
					[clust,exact_text] = getcut([filepath,flnmroot]);
				else
					clust = [];
				end % if strcmpi(button,'yes')
			end % if numSpikesInCutFile ~= numSpikesInTetrode
			idx = strfind(exact_text,': '); % should contain two values
			flist = textscan(exact_text(idx(1)+2:idx(2)-8),'%s','delimiter',',');
			fileStruct = sortstruct(fileStruct,'flnmroot',flist{1});
		else
			clust = [];
		end % if exist([filepath,fileStruct(1).flnmroot,'_',num2str(current_tet),'.cut'],'file')
		mtint.tetrode(itet).cut = clust;
	end
	for i = 1:numel(fileStruct)
		filelist(i) = cellstr(fileStruct(i).flnmroot);
	end % for i = 1:numel(fileStruct)
% use the first file as a header
	mtint.flnmroot = filelist;
	mtint.filepath = filepath;
	mtint.header = getDACQHeader([filepath,filelist{1},'.set'],'set');
	mtint.pos.header = headerCheck(filepath,filelist,'pos');
	mtint = postprocess_DACQ_data( mtint );
else
    mtint = readDACQdata ( filepath, fileStruct.flnmroot );
end % if numel(fileStruct) > 1

% this final check asks the user to select a cut file that corresponds to
% the combined cut file for all the trials, it then reads that file and
% re-orders the fileStruct so that the list of files passed back is
% correctly ordered so that the files are loaded in the same order they
% were cut in Tint
% if numel(fileStruct) > 1
%     
%     
%     [cut_file,filepath] = uigetfile('*.cut','Select the COMBINED .cut file...','MultiSelect','off');
%     [~,exact_text] = getcut([filepath,cut_file]);
%     idx = strfind(exact_text,': '); % should contain two values
%     flist = textscan(exact_text(idx(1)+2:idx(2)-8),'%s','delimiter',',');
%     fileStruct = sortstruct(fileStruct,'flnmroot',flist{1});
% end
