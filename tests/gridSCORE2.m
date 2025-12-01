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
%         gdata.g_score = grid score
%
%%%%%%%% Comments
%   21/08/17 created 
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('method','var') || isempty(method) || all(isnan(method(:)))
    method = 'allen';
end

% preallocate measures
g = NaN; % gscore
gdata = struct; % structure of analysis details
gdata.mid_peak = NaN;
gdata.near_peaks = NaN;
gdata.near_peaks_d = NaN;
gdata.central_ring = NaN;
gdata.g_score = NaN;
gdata.wavelength = NaN;
gdata.radius = NaN;
gdata.central_mask = NaN;
gdata.mean_inner_angle = NaN;
gdata.method = method;

im = single(im);
if all(isnan(im(:)))
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate grid score
switch method
    case {1,'allen'}
        %% Method used by Perez-Escobar et al. (2016) Visual landmarks sharpen grid cell metric and confer context specificity to neurons of the medial entorhinal cortex
        % find blobs
        imb = im>0.1;
        blobs = regionprops(imb,'Centroid','Area','PixelIdxList');
        as = [blobs.Area].';
        blobs = blobs(as>10,:);

        % get distance to image centre
        cents = cell2mat({blobs.Centroid}.');
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
        if length(blobs)==7
            % calculate mean inner angle of hexagon
            sorted_cents = cell2mat({blobs.Centroid}.');
            chull = convhull(sorted_cents(:,1),sorted_cents(:,2));
            sorted_peaks = sorted_cents(chull,:);
            ds_peaks = sqrt(sum((sorted_peaks-cent).^2,2));
            as = NaN(length(ds_peaks)-1,1);
            for i = 1:length(ds_peaks)-1
                l1 = ds_peaks(i);
                l2 = ds_peaks(i+1);
                l3 = sqrt(sum((sorted_peaks(i,:)-sorted_peaks(i+1,:)).^2,2));
                anow = acos((l1^2+l2^2-l3^2)/(2*l1*l2));
                as(i) = anow;
            end
            gdata.mean_inner_angle = rad2deg(circ_mean(as));
        end
        
        % calculate mean distance to closest blobs
        mds = mean(ds(2:end));
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

        % correlate
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
            title(sprintf('g = %.2f, s = %.2f, r = %.2f',g,gdata.wavelength,gdata.radius));   
            caxis([0 nanmax(imcent(:))])
            daspect([1 1 1]);
            axis off
            keyboard
        end

    case {2,'wills'}
        %% Method used by Wills et al. (2012) The abrupt development of adult-like grid cell firing in the medial entorhinal cortex
        % find blobs
        imb = im>0.3;
        blobs = regionprops(imb,'Centroid','Area','PixelIdxList');
        as = [blobs.Area].';
        blobs = blobs(as>3,:);

        % get distance to image centre
        cents = cell2mat({blobs.Centroid}.');
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

        % calculate mean inner angle of hexagon
        sorted_cents = cell2mat({blobs.Centroid}.');
        chull = convhull(sorted_cents(:,1),sorted_cents(:,2));
        sorted_peaks = sorted_cents(chull,:);
        ds_peaks = sqrt(sum((sorted_peaks-cent).^2,2));
        as = NaN(length(ds_peaks)-1,1);
        for i = 1:length(ds_peaks)-1
            l1 = ds_peaks(i);
            l2 = ds_peaks(i+1);
            l3 = sqrt(sum((sorted_peaks(i,:)-sorted_peaks(i+1,:)).^2,2));
            anow = acos((l1^2+l2^2-l3^2)/(2*l1*l2));
            as(i) = anow;
        end
        gdata.mean_inner_angle = rad2deg(circ_mean(as));        
        
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

        % rotate and correlate
        h60 = imrotate(imcent,60,'bilinear','crop');
        h120 = imrotate(imcent,120,'bilinear','crop');
        h30 = imrotate(imcent,30,'bilinear','crop');
        h90 = imrotate(imcent,90,'bilinear','crop');
        h150 = imrotate(imcent,150,'bilinear','crop');

        % correlate
        r1 = corrcoef(h60(:),imcent(:),'rows','pairwise');
        r1 = r1(1,2);
        r2 = corrcoef(h120(:),imcent(:),'rows','pairwise');
        r2 = r2(1,2);
        r3 = corrcoef(h30(:),imcent(:),'rows','pairwise');
        r3 = r3(1,2);
        r4 = corrcoef(h90(:),imcent(:),'rows','pairwise');
        r4 = r4(1,2);
        r5 = corrcoef(h150(:),imcent(:),'rows','pairwise');
        r5 = r5(1,2);
        g = nanmin([r1,r2]) - nanmax([r3,r4,r5]);

        % collect data
        gdata.mid_peak = blobs(1).Centroid;
        gdata.near_peaks = cell2mat({blobs(2:end).Centroid}.');
        gdata.near_peaks_d = ds(2:end);
        gdata.central_ring = imcent;

        gdata.g_score = g;
        gdata.wavelength = mds;
        gdata.radius = sqrt(blobs(1).Area)./pi;

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
            title(sprintf('g = %.2f, s = %.2f, r = %.2f',g,gdata.wavelength,gdata.radius));   
            caxis([0 nanmax(imcent(:))])
            daspect([1 1 1]);
            axis off
            keyboard
        end
end



















        



