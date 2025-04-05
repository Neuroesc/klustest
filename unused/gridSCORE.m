function [g,gdata] = gridSCORE(im)
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
% preallocate measures
g = NaN; % gscore
gdata = struct; % structure of analysis details
gdata.mid_peak = NaN;
gdata.near_peaks = NaN;
gdata.near_peaks_d = NaN;
gdata.central_ring = NaN;
gdata.g_score = NaN;
im = single(im);

if all(isnan(im(:)))
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate grid score
% find blobs
imp = im;
imp(isnan(im)) = 0;
imp(im < 0.2*nanmax(im(:))) = 0;
imat = imregionalmax(imp,26);
[rr,cc] = find(imat);

% distance transform
dmat = zeros(size(im));
dmat(round(size(dmat,1)/2),round(size(dmat,2)/2)) = 1;
dmat = bwdist(dmat); % distance transform of dmat
dd = dmat(sub2ind(size(dmat),rr,cc));
pks = [rr(:),cc(:),dd(:)];

% filter maxima
[~,di] = sort(dd,'ascend');
pks = pks(di,:);
if numel(dd) > 7
    pks = pks(1:7,:);
end
if numel(dd) == 1
    return;
end

% extend peaks
bmat = zeros(size(imat));
bmat(sub2ind(size(imat),pks(:,1),pks(:,2))) = 1;
for pp = 1:numel(pks(:,1))
    val = im(pks(pp,1),pks(pp,2));
    imt = im;
    imt(imt < 0.5*val) = 0;
    imt(imt > 0) = 1;
    ls = bwlabel(imt,8);
    lval = ls(pks(pp,1),pks(pp,2));
    bmat(ls == lval) = 1;
end

% convex hull
hmat = bwconvhull(bmat);
cmat = im;
cmat(~hmat) = NaN;

% get new distance measures to central peak, not central pixel
dmat = zeros(size(im));
dmat(pks(1,1),pks(1,2)) = 1;
dmat = bwdist(dmat); % distance transform of dmat
pks(:,3) = dmat(sub2ind(size(dmat),pks(:,1),pks(:,2)));
[~,di] = sort(pks(:,3),'ascend');
pks = pks(di,:);

% obscure central peak
cmat(dmat < nanmedian(pks(:,3))*0.5) = NaN;
cmat(dmat > nanmedian(pks(:,3))*1.5) = NaN;

% cut to middle, centred on central peak
d2 = double(ceil(nanmedian(pks(:,3))*1.5));
cmatc = padarray(cmat,[d2 d2],NaN,'both');
cmatc = cmatc(pks(1,1)+d2-d2:pks(1,1)+d2+d2,pks(1,2)+d2-d2:pks(1,2)+d2+d2); % cut cmat to a square centred on middle peak

% rotate and correlate
h60 = imrotate(cmatc,60,'bilinear','crop');
h120 = imrotate(cmatc,120,'bilinear','crop');
h30 = imrotate(cmatc,30,'bilinear','crop');
h90 = imrotate(cmatc,90,'bilinear','crop');
h150 = imrotate(cmatc,150,'bilinear','crop');

% correlate
r1 = corrcoef(h60(:),cmatc(:),'rows','pairwise');
r1 = r1(1,2);
r2 = corrcoef(h120(:),cmatc(:),'rows','pairwise');
r2 = r2(1,2);
r3 = corrcoef(h30(:),cmatc(:),'rows','pairwise');
r3 = r3(1,2);
r4 = corrcoef(h90(:),cmatc(:),'rows','pairwise');
r4 = r4(1,2);
r5 = corrcoef(h150(:),cmatc(:),'rows','pairwise');
r5 = r5(1,2);
g = nanmin([r1,r2]) - nanmax([r3,r4,r5]);

% collect data
gdata.mid_peak = pks(1,1:2);
gdata.near_peaks = pks(2:end,1:2);
gdata.near_peaks_d = nanmedian(pks(2:end,3));
gdata.central_ring = cmat;
gdata.g_score = g;

% create figure if required
if 0
    figure
    imc = imagesc(im);
    daspect([1 1 1]);
    msk = ones(size(im)).*0.2;
    msk(~isnan(cmat)) = 1;
    set(imc,'alphadata',msk);
    hold on
    plot(pks(:,2),pks(:,1),'kx','MarkerSize',10);
    title(sprintf('g = %.2f, s = %.2f',g,gdata.near_peaks_d));   
    caxis([0 nanmax(cmat(:))])
    axis off
end
        
        



























        



