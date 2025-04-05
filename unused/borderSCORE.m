function [g,gdata] = borderSCORE(im)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an autocorrelogram image and calculates various grid metrics
%   [g,gdata] = gridSCORE(im)
%
%%%%%%%% Inputs
%   im = the autocorrelogram image
%
%%%%%%%% Outputs
%   g = grid score
%   gdata = structure containing more detailed info
%         gdata.mid_peak = middle peak location [x,y]
%         gdata.near_peaks = surrounding 6 peak locations [x,y]
%         gdata.near_peaks_d = median distance to surrounding peaks
%         gdata.central_ring = autocorrelogram image cut to central fields (with middle removed)
%         gdata.orientation = the angle/orientation of the grid
%         gdata.g_score = grid score
%
%%%%%%%% Comments
%   21/08/17 created 
%   12/09/17 added Allen and Wills methods
%   13/09/17 added peak method
%   15/10/17 removed buggy peak method and corrected Allen method, added orientation
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
% preallocate measures
gdata = struct; % structure of analysis details
gdata.border_score = method;

im = single(im);
if all(isnan(im(:)))
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate grid score       
    % find central peak
    im2 = im;
    im2(im<0.5 | isnan(im) | isinf(im)) = 0; % remove low values and NaNs etc
    cdata = regionprops(logical(im2),'Centroid','EquivDiameter','MajorAxisLength','MinorAxisLength'); % get diameter of blobs

figure
imagesc(im2)
daspect([1 1 1])
keyboard
    
    
    
    
    % get distance to image centre
    cents = cell2mat({cdata.Centroid}.');
    if isempty(cents)
        return
    end        
    cent = [size(im,2)/2,size(im,1)/2];
    ds = sqrt(sum((cents-cent).^2,2)); % all distances to centre of image        
    [~,mindx] = nanmin(ds); % find blob closest to centre

    % get info about this central blob
    crad = cdata(mindx).EquivDiameter/2; % get radius of central blob
    rad_width = crad;
    ccent = cdata(mindx).Centroid; % get centroid of central blob
    min_dist = crad*1.5; % our minimum distance is larger than this
    gdata.border_score = cdata(mindx).MajorAxisLength ./ cdata(mindx).MinorAxisLength; % border score as elongation of central field




















        



