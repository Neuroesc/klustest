function [datout,nfields,thresh_ratemap] = getPFIELDS2(ratemap,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
%getPFIELDS  this function analyses a ratemap for place fields
%
% USAGE:
%           [datout,nfields,thresh_ratemap] = getPFIELDS2(ratemap) process firing rate map with default settings
%
%           [datout,nfields,thresh_ratemap] = getPFIELDS2(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%           of the process
%
%           Parameters include:
%
%           'frcut'           -   (default = 0.2) Scalar, proportion, the firing rate cutoff (proportion of max(ratemap(:))), regions above this may be counted as place fields
%
%           'peakfr'          -   (default = 1) Scalar, Hz, fields with a peak firing rate less than this will be excluded
%
%           'meanfr'          -   (default = 0) Scalar, Hz, fields with a mean firing rate less than this will be excluded
%
%           'minfr'           -   (default = 0) Scalar, Hz, minimum firing rate in Hz - ratemap will be thresholded using frcut and this value
%
%           'arcut'           -   (default = 9) Scalar, pixels, area cutoff the minimum area of contiguous pixels a region must be to count as a place field
%
% OUTPUT:
%           datout - output table, each row corresponds to a place field, see regionprops for a more detailed explanation of each output
%               datout.Area - the area of each field in pixels or ratemap bins
%               datout.Centroid - the centroid of each place field (i.e. the center of the field's shape)
%               datout.WeightedCentroid - the weighted centroid of each field (i.e. it's center of mass when taking into account firing rate)
%               datout.MaxIntensity - the maximum firing rate of each place field
%               datout.MeanIntensity - the mean firing rate of each place field
%               datout.MajorAxisLength - the length of the field's longest axis
%               datout.MinorAxisLength - the length of the field's shortest axis
%               datout.Orientation - the orientation of the field (i.e. the angle between the x-axis and the major axis)
%
%           nfields - number of fields detected
%
%           thresh_ratemap - ratemap thresholded according to input criteria
%
% EXAMPLES:
%     % find the place fields in a ratemap using default settings, plot it, plot the binarized ratemap, add text on top of each field giving its area
%     ratemap = abs(peaks(128)); % fake ratemap
%     [datout,nfields,bmap] = getPFIELDS2(ratemap);
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
% version 6.0.0, Release 02/04/20 changed for quadtest, renamed getPFIELDS2 to avoid confusion, added name value pair arguments so expansion of criteria is possible
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2020 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INPUT ARGUMENTS CHECK
%% Prepare default settings
    def_frcut                   = 0.2;        
    def_minfr                   = 0;             
    def_peakfr                  = 1;    
    def_meanfr                  = 0;        
    def_arcut                   = 9;             
    
%% Parse inputs
    p = inputParser;
    addRequired(p,'ratemap',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    addParameter(p,'frcut',def_frcut,@(x) isnumeric(x) && isscalar(x));  
    addParameter(p,'minfr',def_minfr,@(x) isnumeric(x) && isscalar(x));      
    addParameter(p,'peakfr',def_peakfr,@(x) isnumeric(x) && isscalar(x));  
    addParameter(p,'meanfr',def_meanfr,@(x) isnumeric(x) && isscalar(x));      
    addParameter(p,'arcut',def_arcut,@(x) isnumeric(x) && isscalar(x));  
    parse(p,ratemap,varargin{:});

%% Retrieve parameters 
    config = p.Results;

%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
    % threshold the ratemap using the minimum firing rate proportion (frcut) and the minimum cutoff (minfr)
    thresh_ratemap = imbinarize(ratemap,max([config.minfr, config.frcut*nanmax(ratemap(:))]) );

    % run regionprops to get properties of remaining regions
    datout = regionprops('table',thresh_ratemap,ratemap,'Area','Centroid','WeightedCentroid','MaxIntensity','MeanIntensity','MajorAxisLength','MinorAxisLength','Orientation','ConvexHull');

    % remove fields that are too small or that have a peak firing rate less than peakfr or that have a mean firing rate less than meanfr
    nindx = datout.Area(:) < config.arcut | datout.MaxIntensity(:) < config.peakfr | datout.MeanIntensity(:) < config.meanfr;
    datout(nindx,:) = [];

    % add signal to noise as a measure
    datout.signal_to_noise = datout.MaxIntensity(:,:) ./ nanmean(ratemap(:));

    % find the total number of fields
    nfields = length(datout.Area(:,:));








