function [ctype,ctypen,ctypeb] = getCELLTYPE(sdata,uci,part_now)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an sdata structure, a cell identifier and a part name and finds the cell's putative type
%   ctype = getCELLTYPE(uci,part_now);
%
%%%%%%%% Inputs
%   sdata = sdata structure
%   uci = unique cell identifier
%   pname = part name
%
%%%%%%%% Inputs
%   ctype = text describing the cell type
%   ctypen = number describing the cell type:
%       11 = interneuron
%       12 = fast firing interneuron
%       20 = silent pyramidal
%       21 = medium firing pyramidal
%       22 = fast firing pyramidal
%       213 = place cell (i.e. medium pyramidal with an appended 3)
%       if this number also contains a 4 then the cell is also a grid cell
%       if this number also contains a 5 then the cell is also a HD cell
%
%   ctypeb = binary number describing cell type
%       00000 = silent cell
%       10000 = 0 = interneuron, 1 = pyramidal
%       01000 = 0 = medium frate, 1 = high frate
%       00100 = 0 = low SI, 1 = high SI
%       00010 = 1 = grid cell
%       00001 = 1 = HD cell
%
%%%%%%%% Comments
%   02/04/17 created to contain this code
%   01/11/17 added binary cell typing
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings
place_width_min = 250; % (ms) minimum spoike width allowed for a place cell
place_frate_max = 5; % (Hz) maximum firing rate allowed for a place cell
place_frate_min = 0.1; % (Hz) minimum firing rate allowed for a place cell
place_skagg_min = 0.5; % (bits) minimum spatial information content allowed for a place cell
grid_gscore_min = 0.5; % minimum grid score allowed for a grid cell
hd_rayleigh_min = 0.5; % minimum HD score allowed for a head direction cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine cell type
frate = sdata.(uci).(part_now).frate;
[wpeak,pindx] = nanmax(sdata.(uci).(part_now).waveform_max);
wwide = sdata.(uci).(part_now).waveform_width(pindx);
skagg = sdata.(uci).(part_now).spatial_measures.spatial_information;
gscore = sdata.(uci).(part_now).grid_score;
rvect = sdata.(uci).(part_now).hd_rayleigh;

ctype = [];
ctypen = 00;
ctypeb = NaN;
if wwide > place_width_min % putative pyramidal cells
    if frate > place_frate_max
        ctype = 'fast pyramidal';
        ctypen = 22;
    elseif frate < place_frate_min
        ctype = 'silent cell';
        ctypen = 20;
    else
        if skagg >= place_skagg_min
            ctype = 'place cell';
            ctypen = 213;
        else
            ctype = 'medium pyramidal';
            ctypen = 21;
        end % if skagg > place_skagg_min
    end % if frate > place_frate_cut_max
elseif wwide < place_width_min % putative interneurons
    if frate > place_frate_max
        ctype = 'fast interneuron';
        ctypen = 12;
    else
        ctype = 'interneuron';
        ctypen = 11;        
    end % if frate > place_frate_max
end % if wwide > place_width_min % putative pyramidal cells

if gscore >= grid_gscore_min
    ctype = [ctype ' & grid cell'];
    ctypen = ctypen*10 + 4;
end % gscore >= grid_gscore_min

if rvect >= hd_rayleigh_min
    ctype = [ctype ' & hd cell'];
    ctypen = ctypen*10 + 5;
end % rvect >= hd_rayleigh_min

%% binary method
ctypeb = zeros(1,5);

% 1st input
if wwide > place_width_min
    ctypeb(1) = 1;
end

% 2nd input
if frate < place_frate_min
    return
elseif frate > place_frate_max
    ctypeb(2) = 1;
end

% 3rd input
if skagg >= place_skagg_min
    ctypeb(3) = 1;
end

% 4th input
if gscore >= grid_gscore_min
    ctypeb(4) = 1;
end

% 5th input
if rvect >= hd_rayleigh_min
    ctypeb(5) = 1;
end













