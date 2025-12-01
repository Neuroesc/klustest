function [zmap] = mapDATA3(x,y,z,ls,bs,pr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function
%   [out] = template(in,in2)
%
%%%%%%%% Inputs
%   x = x coordinates of data
%   y = y coordinates of data
%   ls = (default = [min(x),max(x),min(y),max(y)]) the limits you would like in the final map
%       if ls is a single number then ls will equal: min(x)-ls,max(x)+ls,min(y)-ls,max(y)+ls
%   bs = (default = 2) the binsize in cm 
%   pr = (default = 300) the pixel ratio of the data
%
%%%%%%%% Outputs
%   map = the 2d map of the data, unsmoothed
%   xc = a number that can be used to exchange position data back to ratemap coordinates (this would probably need to use the ls numbers as well though)
%
%%%%%%%% Comments
%   08/08/16 created as a simple mapping function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if no limits are given default to the data extremes
if ~exist('ls','var') || isempty(ls)
    x = x + 10;
    y = y + 10;
    ls = [min(x)-10 max(x)+10 min(y)-10 max(y)+10];
elseif numel(ls) == 1 % if the limit is given as a single digit, pad the data using this
    x = x + ls;
    y = y + ls;
    ls = [min(x)-ls max(x)+ls min(y)-ls max(y)+ls];
else % if the limit is given as a vector 
    x = x + ls(1);
    y = y + ls(3);
end % if numel(ls) == 1

% if no binsize is given
if ~exist('bs','var') || isempty(bs)
    bs = 2;
end % if ~exist('bs','var') || isempty(bs)

% if no pixel ratio is given
if ~exist('pr','var') || isempty(pr)
    pr = 300;
end % if ~exist('pr','var') || isempty(pr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate map
bin_size = ceil((bs/100) * pr);
xlength = ceil((ls(2)-ls(1))/bin_size) * bin_size;
ylength = ceil((ls(4)-ls(3))/bin_size) * bin_size;

xbins = xlength/bin_size;
ybins = ylength/bin_size;
nx = xbins;
ny = ybins;

xc = xlength/xbins; % this number can be used to multiply map coordinates and find corresponding position coordinates

%% Do the binning
data = [y x]; 
mapout = hist3(data,{linspace(ls(3),ls(4),ybins) linspace(ls(1),ls(2),xbins)}); % run 2d histogram

[~,~,bindx] = histcounts(x,linspace(ls(1),ls(2),xbins));
[~,~,bindy] = histcounts(y,linspace(ls(3),ls(4),ybins));

count = sparse(bindx,bindy,1,nx,ny);
zsum = sparse(bindx,bindy,z,nx,ny);
k = find(count~=0);
zbar = NaN(size(zsum));
zbar(k) = zsum(k)./count(k);
zmap = zbar;













































