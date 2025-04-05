function [ctype] = getCTYPE(frate,wow,si,g,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A short function to determine the type of neuron, based on a number of characteristics
%
% 25/02/16 created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if frate < 0.1 % If there are few spikes it is noise
    FR = 0;
    ctype = 'noise';
elseif frate > 0.1 && frate < 5 % if the firing rate is between 0.1 and 5 hz it is slow spiking
    FR = 1;
    ctype = 'slow spiking';
elseif frate > 5 % if the firing rate is greater than 5 hz it is fast spiking
    FR = 2;
    ctype = 'fast spiking';
end % if frate < 0.1 % If there are few spikes it is noise

if FR ~= 0
    if wow > 250 % if the width of waveform is greater than 250 it is pyramidal
        ctype = [ctype ', pyramidal'];
    else % otherwise it is an interneuron
        ctype = [ctype ', interneuron'];
    end % if frate < 0.1 % If there are few spikes it is noise

    if si > 0.5 % if the spatial information is greater than 0.5 it is spatially modulated
        ctype = [ctype ', spatial'];
    else % otherwise it is an interneuron
        ctype = [ctype ', non-spatial'];
    end % if frate < 0.1 % If there are few spikes it is noise

    if r > 0.5 % if the rayleigh vector length is greater than 0.5 it is head direction modulated
        ctype = [ctype ', directional'];
    else % otherwise it is an interneuron
        ctype = [ctype ', non-directional'];
    end % if frate < 0.1 % If there are few spikes it is noise

    if g > 0.5 % if the grid score is greater than 0.5 it is griddy
        ctype = [ctype ', grid'];
    else % otherwise it is an interneuron
        ctype = [ctype ', non-grid'];
    end % if frate < 0.1 % If there are few spikes it is noise
end % if FR ~= 0




















end % if FR ~= 0





