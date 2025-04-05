function [btype,ctype,ttype] = getCELLTYPE2(fr,wow,si,gs,rv,fcut,wcut,scut,gcut,rcut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getCELLTYPE2  automated cell identification
% reduces the criteria we use for multiple cell types into one binary output, an integer output and a text description
%
% USAGE:
%         [btype,ctype,ttype] = getCELLTYPE2(fr,wow,si,gs,rv,fcut,wcut,scut,gcut,rcut)
%
% INPUT:
%         fr - cell firing rate
%         wow - cell width of waveform
%         si - spatial information content
%         gs - grid score
%         rv - rayleigh vector length
%         fcut - (optional) firing rate cutoffs [silent,place cell]
%         wcut - (optional) width of waveform cutoff for a place cell
%         scut - (optional) spatial information content cutoff for a place cell
%         gcut - (optional) grid score cutoff for a grid cell
%         rcut - (optional) rayeleigh vector cutoff for a hd cell
%
% OUTPUT:
%    btype - binary output based on the above
%    ctype - cell type as an integer
%    ttype - cell type asd a text description
%
% EXAMPLES:
%
% See also: klustest bin2dec

% HISTORY:
% version 1.0.0, Release 09/04/19 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2019 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
firing_cutoff            = [0.1 5];
wave_width_cutoff        = 250;
spatial_info_cutoff      = 0.5;
grid_score_cutoff        = 0.8;
rayleigh_vector_cutoff   = 0.4;

% deal with input variables
inps = {'fr','wow','si','gs','rv','fcut','wcut','scut','gcut','rcut'};
vals = {'0','0','0','0','0','firing_cutoff','wave_width_cutoff','spatial_info_cutoff','grid_score_cutoff','rayleigh_vector_cutoff'};
reqd = [1 1 1 1 1 0 0 0 0 0];
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
    ttype = 'undefined';
    btype = [true fr>fcut(2) wow>wcut si>scut gs>gcut rv>rcut];
    if fr<fcut(1)
        btype = [true false(1,5)];
    end    
    ctype = bin2dec(num2str(btype,'%d'));

    % start with the most specific and move down to less specific
    if ctype==32
        ttype = 'undefined'; 
    elseif fr<fcut(1)
        ttype = 'silent';    
    elseif btype(5)
        ttype = 'grid cell';
    elseif btype(6)
        ttype = 'HD cell';          
    elseif any(ismember(44:47,ctype))
        ttype = 'place cell';
    elseif any(ismember(60:63,ctype))
        ttype = 'HF place cell';          
    elseif ~btype(3)
        ttype = 'interneuron';  
    elseif btype(3)
        ttype = 'pyramidal';          
    end












































