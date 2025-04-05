function datout = getPFIELDS2(ratemap,frcut,minfr,arcut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getPFIELDS2  this function analyses a ratemap for place fields
% this function takes a ratemap and a thresholding value, thresholds the map and uses regionprops to detect place fields
% and assess their characteristics
%
% USAGE:
%         datout = getPFIELDS2(ratemap,frcut,minfr,arcut)
%
% INPUT:
%         ratemap - a firing rate map like the one produced by mapDATA5
%         frcut - the firing rate cutoff (proportion of max(ratemap(:))), regions above this may be counted as place fields, default is 0.2 or 20%
%         arcut - area cutoff (pixels) the minimum area of contiguous pixels a region must be to count as a place field, default is 9
%
% OUTPUT:
%    datout - output structure
%
% EXAMPLES:
%
% See also: KLUSTEST MAPDATA

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
if ~exist('in','var') || isempty(in)
    in = 
end
stk = dbstack;
function_name = stk.name;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function analyses a ratemap for place fields
%   dataout is a matrix, where each row represents a different field and is composed as follows:
%   [field number, weighted centre x, weighted centre y, centre x, centre y, area (bins), max frate, mean frate, signal to noise ratio, major axis length, minor axis length, orientation, x of max value in field, y of max value in field, max value in field]
%   length(dataout.fields(:,1)) will give the number of fields
%   [datout] = getFIELDS(ratemap)
%
%%%%%%%% Inputs
%   ratemap = the firing rate map to be analysed
%   frcut = Minimum firing rate (% of ratemap max) to be considered a field
%   minfr = (Hz) Minimum cutoff firing rate to be considered a field
%   arcut = The minimum number of pixels to be considered a field, typically people use 9 contiguous pixels
%
%%%%%%%% Outputs
%   dataout is a matrix, where each row represents a different field and is composed as follows:
%   [field number, weighted centre x, weighted centre y, centre x, centre y, area (bins), max frate, mean frate, signal to noise ratio, major axis length, minor axis length, orientation, x of max value in field, y of max value in field, max value in field]
%   length(dataout.fields(:,1)) will give the number of fields
%
%%%%%%%% Comments
%   01/04/15 v1 created
%   22/04/15 v2 created, added regionprops and cleaned up firing thresholds
%   01/04/15 v3 created, fixed problems with field detection
%   08/04/16 v4 created for Dot, added better comments/explanations
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('ratemap','var') || isempty(ratemap)
    error('ERROR: getPFIELDS requires a ratemap... exiting');
end % if ~exist('in','var') || ismepty(in)

if ~exist('frcut','var') || isempty(frcut)
    frcut = 0.2; % Minimum firing rate (% of ratemap max) to be considered a field
end % if ~exist('frcut','var') || ismepty(frcut)

if ~exist('minfr','var') || isempty(minfr)
    minfr = 0; % (Hz) Minimum cutoff firing rate to be considered a field
end % if ~exist('frcut','var') || ismepty(frcut)

if ~exist('arcut','var') || isempty(arcut)
    arcut = 9; % The minimum number of pixels to be considered a field, typically people use 9 contiguous pixels
end % if ~exist('frcut','var') || ismepty(frcut)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise Variables
%fr_cut = nanmean(ratemap(:)); % Minimum firing rate to be considered a field	(based on ratemap mean)
fr_cut = frcut*nanmax(ratemap(:)); % Minimum firing rate to be considered a field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get rid of sub-field firing
ratemap2 = ratemap;	% Duplicate ratemap
ratemap2(ratemap < minfr) = 0;	% make anything less than minimum cutoff firing rate = 0
ratemap2(ratemap < fr_cut) = 0;	% make anything less than minimum cutoff firing rate = 0
ratemap2(isnan(ratemap)) = 0; % make all NaNs = 0 while we are doing this
ratemap2(ratemap2 ~= 0) = 1; % make everything that is not zero = 1

%% Create binary labelled matrix
bin_ratemap = bwlabel(ratemap2,8); % transform ratemap to binary ratemap, each field becomes an island of a single digit representing that field
br2 = bin_ratemap;

%% Localise fields and find their size
stats = regionprops(bin_ratemap,ratemap,'Area','Centroid','WeightedCentroid','MaxIntensity','MeanIntensity','MajorAxisLength','MinorAxisLength','Orientation');% run regionprops to get properties of islands i.e. fields

%% Extract the required data
[pcount,indx] = find([stats.Area] >= arcut); % find which fields have an area over the minimum
br2(~ismember(br2,indx)) = 0;

all_field_data = NaN(length(indx),15); % this will hold the field data temporarily
field_data = NaN(1,15); % preallocate matrix
if ~isempty(indx)
    for i = 1:length(indx) % for every field above min area
        indx1 = indx(i); % get the indice of the current field

        field_data(1:12) = [i stats(indx1).WeightedCentroid stats(indx1).Centroid stats(indx1).Area stats(indx1).MaxIntensity stats(indx1).MeanIntensity (stats(indx1).MaxIntensity/nanmean(ratemap(:))) stats(indx1).MajorAxisLength stats(indx1).MinorAxisLength stats(indx1).Orientation];% extract an align required data in a friendly matrix

        cut_ratemap = ratemap; % duplicate ratemap
        cut_ratemap(bin_ratemap ~= indx1) = 0; % make all values not corresponding to the current field = 0

        [q,p] = find(cut_ratemap == max(nanmax(cut_ratemap))); % find the maximum value in this cut map and it's location
        v(1,1) = ratemap(q(1,1),p(1,1)); % find it's value (assume it's the 1st value in case there are multiples)

        field_data(13:15) = [p(1,1) q(1,1) v(1,1)]; % add this max data to the data matrix
        all_field_data(i,:) = field_data; % accumulate data matrix
    end % for i = 1:length(indx)
end % if ~isempty(indx)
datout.fields = all_field_data; % assign dataout
datout.binary_ratemap = br2;

% %% Here we will make some figures just so you can see what the result is
% figure
% subplot(2,2,1)
% imagesc(ratemap)
% colormap('jet')
% daspect([1 1 1])
% colorbar
% axis xy
% 
% subplot(2,2,2)
% imagesc(bin_ratemap)
% colormap('jet')
% daspect([1 1 1])
% colorbar
% axis xy
% 
% subplot(2,2,3)
% imagesc(ratemap2)
% colormap('jet')
% daspect([1 1 1])
% colorbar
% axis xy
% title(sprintf('%.f fields',length(datout(:,1))))
% 
% for i = 1:length(datout(:,1))
%     hold on
%     plot(datout(i,2),datout(i,3),'g+')
%     text(datout(i,2)+5,datout(i,3)+5,num2str(i),'Color',[1 1 1])
% end % for i = 1:length(datout)
% 
% 








