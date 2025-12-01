function renameforNINGYU(findx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     A stupid function to rename Yave's files and replace some full stops he put in his filenames
%     This function will only replace the two full stops used in the date (the first two full stops in the filename)
%     i.e. 
%     r602_28.06.14_t8_3.cut 
%     would be changed to 
%     r602_28_06_14_t8_3.cut
% 
%     However, the number of full stops replaced can be increased if necessary using the findx input variable
%     i.e. 
%     renameforNINGYU(10)
%     will attempt to replace the first 10 full stops in every file name
% 
%%%%%%%% Inputs
%   findx = (default = 2) how many full stops to replace in every file name (as counted from start)
%
%%%%%%%% Comments
%  07/06/17 created for Ningyu
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('findx','var') || isempty(findx) || all(isnan(findx))
    findx = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_names = dir; % get list of all files in directory
file_names = {file_names.name}.'; % get just their names in a cell array

for nn = 1:length(file_names) % for every file
    curr_file = file_names{nn};
    
    % there are always two dummy filenames, ignore these
    if strcmp(curr_file,'.') || strcmp(curr_file,'..') 
        disp(sprintf('no full stops to replace: %s',curr_file));
        continue
    end
    
    % get an index of all dots
    idx = strfind(curr_file,'.');
    
    % if there is only 1 (for the extension), skip the file
    if numel(idx) == 1 || numel(idx)-1 < findx
        disp(sprintf('no full stops to replace: %s',curr_file));
        continue
    end
    
    % prepare some file names
    new_file = curr_file;
    new_file(idx(1:findx)) = '_';
    file1 = curr_file;					
    file2 = [curr_file, '.temp'];				
    file3 = new_file;			
        
    % actually do the renaming
    disp(sprintf('%s to %s',file1,file3));
    movefile(file1,file2);					
    movefile(file2,file3);					
end




































