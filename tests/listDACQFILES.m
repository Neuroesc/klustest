function [fileStruct,tetsAvailable] = listDACQFILES(filepath,files)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   When a user wants to load multiple files this fcn checks for the
%   presence/ absence of tetrode, pos files etc and re-orders files based on
%   the combined cut file so the order of files loaded is the same as they
%   were cut in Tint.
%   [fileStruct,tetsAvailable] = listDACQFILES()
%
%%%%%%%% Inputs
%   filepath = the path to the files (usually the working directory)
%   files = the file names that were analysed, if multiple sessions were merged then their names must be given as a cell array
%
%%%%%%%% Outputs
%   fileStruct = a structure containing the file info
%   tetsAvailable = the tetrodes that can be analysed
%
%%%%%%%% Comments
%   05/08/16 created from list_all_DACQ_files and 'Rodified'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('filepath','var') || isempty(filepath) || ~exist('files','var') || isempty(files)
    error('No files given to listDACQFILES... exiting')
end % ~exist('filepath','var') || isempty(filepath) || ~exist('files','var') || isempty(files)
if ~strcmp(filepath(end),filesep)
    filepath(end+1) = filesep;
end % if ~strcmp(filepath(end),filesep)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a structure to hold the names of files etc
fileStruct = struct('flnmroot',{},'posFilePresent',{},'tetrodeFiles',{},'eegFiles',{});

for ifile = 1:numel(files)
    % check all pos files are present
    fileStruct(ifile).flnmroot = files{ifile};
    if isempty(dir([filepath,fileStruct(ifile).flnmroot,'.pos']))
        error('filePresenceCheck:posFile',[fileStruct(ifile).flnmroot, '.pos is missing so exiting']);
    else
        fileStruct(ifile).posFilePresent = 1;
    end % if isempty(dir([filepath,fileStruct(ifile).flnmroot,'.pos']))
    
    % check for tetrode files
    tet_filelist = dir([filepath,fileStruct(ifile).flnmroot,'.*']);
    f_type = zeros(numel(files),numel(tet_filelist));
    for i = 1:numel(tet_filelist)
        f_type(i) = str2double(tet_filelist(i).name(strfind(tet_filelist(i).name,'.')+1:end));
    end % for i = 1:numel(tet_filelist)
    f_type(isnan(f_type) | (f_type == 0)) = [];
    fileStruct(ifile).tetrodeFiles = f_type;
    
    % check for eeg files
    if isempty(dir([filepath,fileStruct(ifile).flnmroot,'.eeg*']))
        warning('filePresenceCheck:eegFile',[fileStruct(ifile).flnmroot, '.eeg is missing'])
    else
        fileStruct(ifile).eegFiles = 1;
    end % if isempty(dir([filepath,fileStruct(ifile).flnmroot,'.eeg*']))
end % for ifile = 1:numel(files)

% check what tetrode files are present across trials
tetsAvailable = fileStruct(1).tetrodeFiles;
for ifile = 2:numel(fileStruct)
    tetsAvailable = intersect(tetsAvailable,fileStruct(ifile).tetrodeFiles);
end % for ifile = 2:numel(fileStruct)

