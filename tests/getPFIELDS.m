function [datout,nfields,thresh_ratemap] = getPFIELDS(ratemap,frcut,minfr,arcut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getPFIELDS  this function analyses a ratemap for place fields
% this function takes a ratemap and a thresholding value, thresholds the map and uses regionprops to detect place fields
% and assess their characteristics
%
% USAGE:
%         [datout,nfields,thresh_ratemap] = getPFIELDS(ratemap,frcut,minfr,arcut)
%
% INPUT:
%         ratemap - a firing rate map like the one produced by mapDATA
%         frcut - the firing rate cutoff (proportion of max(ratemap(:))), regions above this may be counted as place fields, default is 0.2 or 20%
%         minfr - minimum firing rate in Hz (i.e. pixels must have a value > frcut * max(ratemap(:)) && minfr), default is 0
%         arcut - area cutoff (pixels) the minimum area of contiguous pixels a region must be to count as a place field, default is 9
%
% OUTPUT:
%    datout - output table, each row corresponds to a place field, see regionprops for a more detailed explanation of each output
%       datout.Area - the are of each field in pixels or ratemap bins
%       datout.Centroid - the centroid of each place field (i.e. the center of the field's shape)
%       datout.WeightedCentroid - the weighted centroid of each field (i.e. it's center of mass when taking into account firing rate)
%       datout.MaxIntensity - the maximum firing rate of each place field
%       datout.MeanIntensity - the mean firing rate of each place field
%       datout.MajorAxisLength - the length of the field's longest axis
%       datout.MinorAxisLength - the length of the field's shortest axis
%       datout.Orientation - the orientation of the field (i.e. the angle between the x-axis and the major axis)
%
% EXAMPLES:
%     % find the place fields in a ratemap using default settings, plot it, plot the binarized ratemap, add text on top of each field giving its area
%     ratemap = abs(peaks(128));
%     [datout,nfields,bmap] = getPFIELDS(ratemap);
%     figure
%     subplot(1,2,1)
%     imagesc(ratemap);
%     daspect([1 1 1])
%     subplot(1,2,2)
%     imagesc(bmap);
%     daspect([1 1 1])
%     title(sprintf('fields: %d',nfields))
%     hold on
%     text(datout.Centroid(:,1),datout.Centroid(:,2),cellstr(num2str(datout.Area(:))),'Color','w')
%
% See also: KLUSTEST REGIONPROPS IMBINARIZE MAPDATA

% HISTORY:
% version 1.0.0, Release 01/04/15 Initial release
% version 2.0.0, Release 22/04/15 added regionprops and cleaned up firing thresholds
% version 3.0.0, Release 01/04/15 fixed problems with field detection
% version 4.0.0, Release 08/04/16 updated, cleaned, added better comments/explanations
% version 5.0.0, Release 08/08/18 overhauled for klustest update, simplified, removed loop, changed to a table output
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with missing inputs
inps = {'frcut','minfr','arcut'};
vals = {'0.2','0','9'};
for ff = 1:length(inps)
    if ~exist(inps{ff},'var')
        eval([inps{ff} '=' vals{ff} ';']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
% threshold the ratemap using the minimum firing rate proportion, frcut, and the minimum cutoff, minfr
thresh_ratemap = imbinarize(ratemap,max([minfr, frcut*nanmax(ratemap(:))]) );

% run regionprops to get properties of remaining regions
datout = regionprops('table',thresh_ratemap,ratemap,'Area','Centroid','WeightedCentroid','MaxIntensity','MeanIntensity','MajorAxisLength','MinorAxisLength','Orientation','ConvexHull');

% remove fields that are too small
nindx = datout.Area(:) < arcut;
datout(nindx,:) = [];

% add signal to noise as a measure
datout.signal_to_noise = datout.MaxIntensity(:,:) ./ nanmean(ratemap(:));

% find the total number of fields
nfields = length(datout.Area(:,:));








