function [g,gdata] = gridSCORE(im,method)
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
if ~exist('method','var') || isempty(method) || all(isnan(method(:)))
    method = 'langston';
end

% preallocate measures
g = NaN; % gscore
gdata = struct; % structure of analysis details
gdata.mid_peak = NaN;
gdata.near_peaks = NaN;
gdata.near_peaks_d = NaN;
gdata.central_ring = NaN;
gdata.g_score = NaN;
gdata.c_score = NaN;
gdata.border_score = NaN;
gdata.g_score2 = NaN;
gdata.wavelength = NaN;
gdata.radius = NaN;
gdata.central_mask = NaN;
gdata.ring_mask = NaN;
gdata.mean_inner_angle = NaN;
gdata.orientation = NaN;
gdata.spokes = NaN;
gdata.method = method;

im = single(im);
if all(isnan(im(:)))
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate grid score
switch method       
    case {'allen'}
        %% Method used by Wills et al. (2012) The abrupt development of adult-like grid cell firing in the medial entorhinal cortex
        % find blobs
        imb = im>0.1;
        blobs = regionprops(imb,'Centroid','Area','PixelIdxList');
        as = [blobs.Area].';
        blobs = blobs(as>10,:);

        % get distance to image centre
        cents = cell2mat({blobs.Centroid}.');
        if isempty(cents)
            return
        end
        cent = [size(im,2)/2,size(im,1)/2];
        ds = sqrt(sum((cents-cent).^2,2)); % all distances to centre of image

        % find middle peak and recalculate distance to it instead
        [~,pindx] = min(ds); % find peak closest to centre - this is the origin
        cent = cents(pindx,:);
        ds = sqrt(sum((cents-cent).^2,2)); % all distances to centre of image

        % sort blobs according to distance
        [ds,sindx] = sort(ds,'ascend');
        blobs = blobs(sindx,:);
        if length(blobs)==1
            return
        end
        if length(blobs)>7
            blobs = blobs(1:7,:);
            ds = ds(1:7,:);
        end       
        
        % calculate grid orientation
        L = createLine(ones(length(blobs)-1,2).*blobs(1).Centroid,cell2mat({blobs(2:end).Centroid}'));
        as = 360 - rad2deg(lineAngle(L));
        grid_ori = min(as);        
        
        % calculate mean distance to closest blobs
        mds = mean(ds);
        dcut = ceil(mds*1.25);
        dcuti = ceil(mds*0.4);

        % cut to the central portion of the autocorrelation
        rcent = round([blobs(1).Centroid(2),blobs(1).Centroid(1)]);
        imp = padarray(im,[dcut dcut],NaN,'both'); % pad array - sometimes mds is calculated diagonally and is larger than im is wide
        imcent = imp(rcent(1)+dcut-dcut:rcent(1)+dcut+dcut,rcent(2)+dcut-dcut:rcent(2)+dcut+dcut); % take the central part of the padded image
        dmat = zeros(size(imcent));
        dmat(ceil(size(imcent,2)/2),ceil(size(imcent,1)/2)) = 1;
        dmat = bwdist(dmat);
        imcent(dmat>dcut) = NaN;

        % remove the central peak
        imcent(dmat<dcuti) = NaN;

        % rotational correlation
        rs = NaN(5,1);
        as = 30:30:150;
        for a = 1:length(as)
            mrot = imrotate(imcent,as(a),'bilinear','crop');
            r = corrcoef(mrot(:),imcent(:),'rows','pairwise'); 
            rs(a) = r(1,2);
        end
        g = ((rs(2)+rs(4))/2) - ((rs(1)+rs(3)+rs(5))/3);

        % collect data
        gdata.mid_peak = blobs(1).Centroid;
        gdata.near_peaks = cell2mat({blobs(2:end).Centroid}.');
        gdata.near_peaks_d = ds(2:end);
        gdata.central_ring = imcent;

        gdata.g_score = g;
        gdata.wavelength = mds;
        gdata.radius = sqrt(blobs(1).Area)./pi;
        gdata.orientation = grid_ori;
        gdata.spokes = L;

        dmat = zeros(size(im));
        dmat(rcent(1),rcent(2)) = 1;
        dmat = bwdist(dmat);
        msk = ones(size(im)).*0.2;
        msk(dmat<dcut & dmat>dcuti) = 1;
        gdata.central_mask = msk;

        % create figure if required
        if 0
            figure
            imc = imagesc(im);
            set(imc,'alphadata',msk);    
            hold on
            plot(gdata.near_peaks(:,1),gdata.near_peaks(:,2),'kx','MarkerSize',10);
            title(sprintf('g = %.2f, s = %.2f, r = %.2f, o = %.2f',g,gdata.wavelength,gdata.radius,gdata.orientation));   
            drawLine(L);
            caxis([0 nanmax(imcent(:))])
            daspect([1 1 1]);
            axis off
            keyboard
        end        
        
    case {'wills'}
        %% Method used by Wills et al. (2012) The abrupt development of adult-like grid cell firing in the medial entorhinal cortex
        % find blobs
        imb = im>0.3;
        blobs = regionprops(imb,'Centroid','Area','PixelIdxList');
        as = [blobs.Area].';
        blobs = blobs(as>3,:);

        % get distance to image centre
        cents = cell2mat({blobs.Centroid}.');
        if isempty(cents)
            return
        end        
        cent = [size(im,2)/2,size(im,1)/2];
        ds = sqrt(sum((cents-cent).^2,2)); % all distances to centre of image

        % find middle peak and recalculate distance to it instead
        [~,pindx] = min(ds); % find peak closest to centre - this is the origin
        cent = cents(pindx,:);
        ds = sqrt(sum((cents-cent).^2,2)); % all distances to centre of image

        % sort blobs according to distance
        [ds,sindx] = sort(ds,'ascend');
        blobs = blobs(sindx,:);
        if length(blobs)==1
            return
        end
        if length(blobs)>7
            blobs = blobs(1:7,:);
            ds = ds(1:7,:);
        end       
        
        % calculate grid orientation
        L = createLine(ones(length(blobs)-1,2).*blobs(1).Centroid,cell2mat({blobs(2:end).Centroid}'));
        as = 360 - rad2deg(lineAngle(L));
        grid_ori = min(as);        
        
        % calculate mean distance to closest blobs
        mds = mean(ds);
        dcut = ceil(mds*1.25);
        dcuti = ceil(mds*0.4);

        % cut to the central portion of the autocorrelation
        rcent = round([blobs(1).Centroid(2),blobs(1).Centroid(1)]);
        imp = padarray(im,[dcut dcut],NaN,'both'); % pad array - sometimes mds is calculated diagonally and is larger than im is wide
        imcent = imp(rcent(1)+dcut-dcut:rcent(1)+dcut+dcut,rcent(2)+dcut-dcut:rcent(2)+dcut+dcut); % take the central part of the padded image
        dmat = zeros(size(imcent));
        dmat(ceil(size(imcent,2)/2),ceil(size(imcent,1)/2)) = 1;
        dmat = bwdist(dmat);
        imcent(dmat>dcut) = NaN;

        % remove the central peak
        imcent(dmat<dcuti) = NaN;

        % rotational correlation
        rs = NaN(5,1);
        as = 30:30:150;
        for a = 1:length(as)
            mrot = imrotate(imcent,as(a),'bilinear','crop');
            r = corrcoef(mrot(:),imcent(:),'rows','pairwise'); 
            rs(a) = r(1,2);
        end
        g = nanmin([rs(2),rs(4)]) - nanmax([rs(1),rs(3),rs(5)]);        

        % collect data
        gdata.mid_peak = blobs(1).Centroid;
        gdata.near_peaks = cell2mat({blobs(2:end).Centroid}.');
        gdata.near_peaks_d = ds(2:end);
        gdata.central_ring = imcent;

        gdata.g_score = g;
        gdata.wavelength = mds;
        gdata.radius = sqrt(blobs(1).Area)./pi;
        gdata.orientation = grid_ori;
        gdata.spokes = L;

        dmat = zeros(size(im));
        dmat(rcent(1),rcent(2)) = 1;
        dmat = bwdist(dmat);
        msk = ones(size(im)).*0.2;
        msk(dmat<dcut & dmat>dcuti) = 1;
        gdata.central_mask = msk;

        % create figure if required
        if 0
            figure
            imc = imagesc(im);
            set(imc,'alphadata',msk);    
            hold on
            plot(gdata.near_peaks(:,1),gdata.near_peaks(:,2),'kx','MarkerSize',10);
            title(sprintf('g = %.2f, s = %.2f, r = %.2f, o = %.2f',g,gdata.wavelength,gdata.radius,gdata.orientation));   
            drawLine(L);
            caxis([0 nanmax(imcent(:))])
            daspect([1 1 1]);
            axis off
            keyboard
        end
        
    case {'langston'}
        %% Method based on that used by Langston et al. (2010) Development of the Spatial Representation System in the Rat
        % rings are cut from the autocorrelation at different distances from the centre and the standard rotation and correlation
        % is performed on each one. The highest of these grid scores is used, the radius of the ring is the grid spacing. A sine wave
        % is fitted to the values in the ring and this is used to estimate the positions of the fields. The grid orientation is the
        % angle from the centre to the first one of these fields, counter-clockwise. The grid field radius is taken as the radius of
        % the centre field.
        rad_steps = 1; % amount by which to increase radii each step
        
        % find central peak
        im2 = im;
        im2(im<0.3 | isnan(im) | isinf(im)) = 0; % remove low values and NaNs etc
        cdata = regionprops(logical(im2),'Centroid','EquivDiameter','MajorAxisLength','MinorAxisLength'); % get diameter of blobs
        
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
        
        % distance matrix
        dmat = zeros(size(im)); % prepare a matrix of zeros, same size as autocorr
        dmat(ceil(ccent(2)),ceil(ccent(1))) = 1; % centre blob equals 1        
        dmat = bwdist(dmat); % every pixel is the distance to this
        
        %% calculate grid score on ever increasing radii
        ds = min_dist:rad_steps:((max(size(im))/2)-crad); % distances we will use
        if isempty(ds) % if the central peak is too large we can't get a ring anywhere
            return
        end
        
        gs = NaN(length(ds),1); % prepare an empty vector
        for dd = 1:length(ds) % for every distance we want to test
            d2 = dmat; % get a copy of the distance matrix
            i2 = im; % get a copy of the autocorr
            i2(d2<ds(dd)-rad_width | d2>ds(dd)+rad_width) = NaN; % cut to a doughnut of radius dd and width rad_width*2
            
            % rotational correlation
            rs = NaN(5,1); % prepare empty vector
            as = 30:30:150; % angles we want to test
            for a = 1:length(as) % for every angle
                mrot = imrotate(i2,as(a),'bilinear','crop'); % rotate the autocorr by this angle
                r = corrcoef(mrot(:),i2(:),'rows','pairwise'); % correlate the rotated map with the original
                rs(a) = r(1,2); % collect the r value
            end
            gs(dd) = nanmin([rs(2),rs(4)]) - nanmax([rs(1),rs(3),rs(5)]); % gs is the grid score    
        end
        [g,dindx] = nanmax(gs); % get the maximum grid score and find the distance associated with it
        
        %% calculate the orientation of the grid
        pref_width = ds(dindx); % pref_width is the distance with the best grid score
        gspacing = pref_width; % the spacing of the grid fields equals the radius of this doughnut
        msk = true(size(im)); % prepare a mask for cutting to this doughnut
        msk(dmat<pref_width-rad_width/2 | dmat>pref_width+rad_width/2) = false; % remove everything outside the doughnut
        im3 = im; % get a copy of the autocorr
        im3(~msk) = NaN; % cut it to the doughnut with the best grid score   

        [cc,rr] = meshgrid(1:size(im,2),1:size(im,1)); % indices of every pixel in autocorr
        as = atan2d(rr-ceil(ccent(2)),cc-ceil(ccent(1))); % convert these to angles to the centre blob
        as(dmat<pref_width-rad_width/2 | dmat>pref_width+rad_width/2) = NaN; % cut to doughnut
        as = abs(wrapTo360(as)-360); % convert to degrees
        
        % fit function to data in order to find peaks and determine orientation
        rindx = isnan(as(:)) | isnan(im3(:)) | isinf(as(:)) | isinf(im3(:));
        asb = double(as(:));
        asb(rindx) = [];
        im3b = double(im3(:));
        im3b(rindx) = [];
            
        % smoothing spline fit
        if numel(asb) < 10
            g2 = NaN;
            grid_ori = NaN;
        else
            [xt,yt] = prepareCurveData(asb,im3b);
            yt = yt-mean(yt(:));
            ft = fittype('smoothingspline');
            opts = fitoptions('Method','SmoothingSpline');
            opts.SmoothingParam = 0.04718;
            [fitresult,gof] = fit(xt,yt,ft,opts);

            % goodness of fit and curve values
            g2 = gof.rsquare;
            xn = 0:1:360;
            yn = fitresult(xn);

            % find peaks
            [~,grid_ori] = findpeaks(yn,xn,'NPeaks',1); 
            if isempty(grid_ori)
                grid_ori = NaN;
            end
        end
        
        % accumulate data
        gdata.orientation = grid_ori;
        gdata.g_score2 = g2;
        gdata.mid_peak = ccent;
        gdata.g_score = g;
        gdata.wavelength = gspacing;
        gdata.radius = crad;
        gdata.ring_mask = msk;

        % create figure if required
        if 0
            figure
            imc = imagesc(im);
            msk(~msk) = 0.2;
            set(imc,'alphadata',msk);    
            hold on
            title(sprintf('g = %.2f, s = %.2f, r = %.2f, o = %.2f',g,gdata.wavelength,gdata.radius,gdata.orientation));   
            caxis([0 nanmax(im(:))])
            daspect([1 1 1]);
            axis off
            keyboard
        end        
        
    case {'soman'} 
        % Soman, K., Muralidharan, V., and Chakravarthy, V.S. (2017)
        % A Model of Multisensory Integration and its Influence on Hippocampal Spatial Cell Responses
        % IEEE Transactions on Cognitive and Developmental Systems
        
        % rotational correlation
        rs = NaN(5,1);
        as = 30:30:150;
        for a = 1:length(as)
            mrot = imrotate(im,as(a),'bilinear','crop');
            r = corrcoef(im(:),mrot(:),'rows','pairwise'); 
            rs(a) = r(1,2);
        end
        g = nanmin([rs(2),rs(4)]) - nanmax([rs(1),rs(3),rs(5)]);           
        
        % accumulate data
        gdata.g_score = g;        
        
    case {'grieves'} 
        % threshold autocorrelation and detect fields
        abin = imbinarize(im,0.3);
        if ~any(abin)
            return
        end
        sts = regionprops('table',abin,'centroid','EquivDiameter','MajorAxisLength');
        sts.ds = sqrt(sum((sts.Centroid-size(abin,[2 1])./2).^2,2));
        sts = sortrows(sts,4);
        central_peak = sts.Centroid(1,:);
        
        % find extent of central field
        abin = imbinarize(im,0.1);
        sts = regionprops('table',abin,'centroid','EquivDiameter','MajorAxisLength');   
        sts.ds = sqrt(sum((sts.Centroid-size(abin,[2 1])./2).^2,2));
        sts = sortrows(sts,4);        
        gdata.border_score = sts.MajorAxisLength(1,:) ./ hypot(size(im,1)/2,size(im,2)/2); % length of the central blob as proportion of autocorr width
        
        central_diameter = sts.EquivDiameter(1,:);
        if size(sts,1)<7
            inner_circle = sts.Centroid(2:end,:);            
        else
            inner_circle = sts.Centroid(2:7,:);
        end
        
        % get their angles and mean distance
        [th,ra] = cart2pol( inner_circle(:,1)-central_peak(1), inner_circle(:,2)-central_peak(2) );
        [inner_circle_pol,i] = sortrows([th ra],2);
        inner_circle = inner_circle(i,:);
        
%         % correct the larger radii to match the average of the others
%         inner_circle_pol(end-1:end,2) = nanmean(inner_circle_pol(1:end-2,2));
%         [x,y] = pol2cart(inner_circle_pol(:,1),inner_circle_pol(:,2));
%         inner_circle2 = [x,y] + central_peak;
%         
%         % apply the affine transformation as long as the scaling and
%         % shearing factors are not too large
%         tf = fitgeotrans([inner_circle;central_peak],[inner_circle2;central_peak],'affine');
%         if ~any( abs( [1 0 0 1]-tf.T([1 2 4 5]) ) > 0.15 )
%             im2 = imwarp(im,tf,'OutputView',imref2d(size(im)));
%         else
%             inner_circle2 = inner_circle;
%             inner_circle_pol = [th ra];
            im2 = im;
%         end

        % make a mask removing the central peak and cutting to the outer 6
        % peaks
        ik = zeros(size(im2));
        ik(round(central_peak(2)),round(central_peak(1))) = 1;
        ib = bwdist(ik);
        msk = ib < nanmean(inner_circle_pol(:,2))+central_diameter/2 & ib > nanmean(inner_circle_pol(:,2))-central_diameter/2 & ib > central_diameter/2;
        im3 = im2;
        im3(~msk) = NaN;
        
        % get inner angles of fields
        angs = sort(th(:));
        a1 = angdiff(angs(:), circshift(angs(:),1));
        a1_d60 = abs( a1 - deg2rad(60) );
        s60 = circ_std(a1_d60);
        m60 = circ_mean(a1_d60);
        ori = rad2deg( min(angs(angs>0)) );

        % rotational correlation
        rs = NaN(5,1);
        as = 30:30:150;
        for a = 1:length(as)
            mrot = imrotate(im3,as(a),'nearest','crop');
            rs(a) = corr(im3(:),mrot(:),'type','Pearson','rows','pairwise'); 
        end
        g = nanmin([rs(2),rs(4)]) - nanmax([rs(1),rs(3),rs(5)]);  
        
        % rotational correlation
        rs = NaN(5,1);
        as = 0:45:180;
        for a = 1:length(as)
            mrot = imrotate(im3,as(a),'nearest','crop');
            rs(a) = corr(im3(:),mrot(:),'type','Pearson','rows','pairwise'); 
        end
        c = nanmin([rs(1),rs(3),rs(5)]) - nanmax([rs(2),rs(4)]);         
        
        % accumulate data
        gdata.mid_peak = central_peak;
        gdata.near_peaks = inner_circle;
        gdata.near_peaks_d = inner_circle_pol(:,2);
        gdata.central_ring = nanmean(inner_circle_pol(:,2));
        gdata.g_score = g;
        gdata.c_score = c;   
        gdata.angle_deviation_around_60 = s60;
        gdata.angle_mean_diff_60 = m60;
        gdata.wavelength = nanmean(inner_circle_pol(:,2));
        gdata.radius = central_diameter/2;
        gdata.central_mask = ib < central_diameter/2;
        gdata.ring_mask = msk;
%         gdata.mean_inner_angle = adiff;
        gdata.orientation = ori;
        gdata.method = method;
% 
% figure
% subplot(2,3,1)
% imagesc(im); 
% daspect([1 1 1])
% 
% subplot(2,3,2)
% imagesc(abin); 
% daspect([1 1 1])
% 
% subplot(2,3,3)
% imagesc(abin); hold on;
% daspect([1 1 1])
% plot(inner_circle(:,1),inner_circle(:,2),'wx')
% plot(inner_circle2(:,1),inner_circle2(:,2),'gx')
% plot(central_peak(:,1),central_peak(:,2),'rx')
% 
% subplot(2,3,4)
% imagesc(im2); 
% daspect([1 1 1])
% 
% subplot(2,3,5)
% img = imagesc(im2); 
% set(img,'alphadata',(double(msk)+1)/2);
% daspect([1 1 1])
% keyboard
  
end



















        



