function [strains,offsets] = shiftTRAIN(spt,ts,nshuff,offlims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function
%   [out] = template(in,in2)
%
%%%%%%%% Inputs
%   spt = spike times (s)
%   ts = (default = [min(spt) max(spt)]) vector of time limits [min max], this should be session min and max (s)
%   nshuff = (default = 1) the number of spike shuffles to do
%   offlims = (default = [-20 20]) the minimum and maximum range of time to shift spikes (s)
%
%%%%%%%% Outputs
%   strains = shifted spike trains, one spike per row, one shift per column
%   offsets = the shifts used to shift spike trains
%
%%%%%%%% Comments
%   17/05/17 created 
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
% if no number of shuffles was provided
if ~exist('nshuff','var') || isempty(nshuff) || all(isnan(nshuff))
    nshuff = 1;
end

% if no session time limits were provided
if ~exist('ts','var') || isempty(ts) || all(isnan(ts))
    ts = [min(spt) max(spt)];
    disp(sprintf('WARNING: shiftTRAIN using default session limits: %.1f - %.1f...',ts(1),ts(2)));
end

% if no offset limits were provided
if ~exist('offlims','var') || isempty(offlims) || all(isnan(offlims))
    offlims = [-20 20];
end
offsets = offlims(1) + (offlims(2)-offlims(1)) .* rand(nshuff,1);

spt = double(spt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate shuffle(s)
strains = NaN(numel(spt),nshuff);
for ii = 1:nshuff
    % normalise so minimum spt is zero (x is now distance from mn)
    vi2 = spt - ts(1);
    mn2 = 0;
    mx2 = ts(2) - ts(1);

    % enforce offset
    vi3 = vi2 + offsets(ii,1);

    % remove multiple offsets (shifting by 3 lengths is functionally the same as moving by one)
    vi3 = vi3 - fix(min(vi3) / mx2) * mx2;

    % rotate extra maximum values to beginning and extra minimum values to end
    vi3(vi3 > mx2) = vi3(vi3 > mx2) - mx2;
    vi3(vi3 < mn2) = mx2 - vi3(vi3 < mn2);

    % return vector from normalisation
    vi4 = vi3 + ts(1);
    vi4 = sort(vi4,'ascend');

    % collect data
    strains(:,ii) = vi4(:);
end



























