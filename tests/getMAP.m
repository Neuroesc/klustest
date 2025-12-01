function [mapout,xchange] = getMAP(x,y,maxx,maxy,pr,bs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a simple function which takes x and y data and bins it
% bs = the size of the bins in cm
% pr = the number of pixels per metre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup initial variables
bin_size = ceil((bs/100) * pr);
xlength = ceil(maxx/bin_size)*bin_size;
ylength = ceil(maxy/bin_size)*bin_size;

xbins = xlength/bin_size;
ybins = ylength/bin_size;

xchange = xlength/xbins;														% this number can be used to multiply map coordinates and find corresponding position coordinates

%% Do the binning
data = [y x]; 
mapout = hist3(data,{linspace(0,ylength,ybins) linspace(0,xlength,xbins)});                                                          % run 2d histogram
nindx = isnan(mapout);
mapout(nindx) = 0;



































