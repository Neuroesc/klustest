function pdata = phasePRECESS2D(sdata,uci,part_now)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an sdata structure from Klustest and performs the 2D phase precession analysis outlined by
%   Jeewajee et al. (2013) Theta phase precession of grid and place cell firing in open environments
%   Field boundaries are detected, trajectories through these fields are found, these are then rotated and projected onto
%   a horizontal line passing left-right. Spike phase can then be expressed in termns of position through the field
%   We then perform a non-parametric linear-circular regression on this data to find a relationship between the two (code from Zugaro)
%   Obviously if there are no distinct place or grid fields this function will not perform well.
%
%   pdata = phasePRECESS2D(sdata,uci,part_now)
%
%%%%%%%% Inputs
%   sdata = sdata structure
%   uci = unique cell identifier
%   part_now = part name
%
%%%%%%%% Outputs
%   pdata = structure containing phase precession data
%         pdata.n_trajectories = a vector, one row per place field, giving the number of trajectories through that field
%         pdata.regression = vector, regression slope, intercept and coefficient
%         pdata.theil_regression = vector, regression slope, intercept and coefficient for Theil–Sen regression
%         pdata.regression_points = matrix, first column is the position of every spike relative to its field, second column is the spike phase
%
%%%%%%%% Comments
%   01/10/17 created to test 2D phase precession, hope to extend to 3D
%   14/11/17 renamed phasePRECESS2D and added to klustest
%   16/11/17 replaced Zugaro's circular-linear regression with my own
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
min_fsize = sdata.config.arcut; % min size for place fields
fcut = sdata.config.frcut; % percentage cutoff firing rate for detecting fields

pdata = struct;
pdata.n_trajectories = NaN;
    pdata.regression_info.slope = NaN;
    pdata.regression_info.phase = NaN;
    pdata.regression_info.slope_deg = NaN;
    pdata.regression_info.phase_deg = NaN;
    pdata.regression_info.fit_rmse = NaN;
    pdata.regression_info.fit_mean = NaN;
    pdata.regression_info.circ_corr = [NaN,NaN];
    pdata.regression_info.img = zeros(128,128);
    pdata.regression_info.img_slope = NaN;
    pdata.regression_info.img_phase = NaN;
    pdata.regression_info.img_phase_offset = NaN;
    pdata.regression_info.img_res = 128;
pdata.regression_points = [NaN,NaN];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load ratemap and position data
ratemap = sdata.(uci).(part_now).ratemap;
pox = sdata.(uci).(part_now).ratemap_data.poxnew;
poy = sdata.(uci).(part_now).ratemap_data.poynew;
pot = sdata.(part_now).pot;
pov = sdata.(part_now).pov;
spx = sdata.(uci).(part_now).ratemap_data.spxnew;
spy = sdata.(uci).(part_now).ratemap_data.spynew;
spt = sdata.(uci).(part_now).spt;
spp = sdata.(uci).(part_now).spike_phase;

%% detect fields
rmap = ratemap;
rmap(rmap<fcut*nanmax(rmap(:)) | isnan(rmap)) = 0;
rmap(rmap>0) = 1;
sts = regionprops(logical(rmap),ratemap,'Area','Centroid','MaxIntensity','ConvexHull','PixelIdxList','EquivDiameter');

%% For every field
% disp(sprintf('%d fields',length(sts)));
pmat = [];
smat = [];
ntraj = NaN(length(sts),1);
for ff = 1:length(sts) % for every place field
    % get field details
    farea = sts(ff).Area;
    if farea < min_fsize
        continue
    end
    chull = sts(ff).ConvexHull;
    pindx = sts(ff).PixelIdxList;
    ctroi = sts(ff).Centroid;
    ediam = sts(ff).EquivDiameter;
    
    % cut map to field
    mapnow = zeros(size(ratemap));
    mapnow(pindx) = ratemap(pindx);
    
    % find trajectories passing through it
    tindx = inpolygon(pox,poy,chull(:,1),chull(:,2));
    tindx = bwlabel(tindx);
    %disp(sprintf('%d trajectories through field',max(tindx)));
    ntraj(ff,1) = max(tindx);

    %% For every trajectory
    for tt = 1:max(tindx)
        if sum(tindx==tt) < 20
            continue
        end
        % current trajectory
        pnow = [pox(tindx==tt),poy(tindx==tt)];
        vnow = pov(tindx==tt)*100; % instantaneous velocity in cm/s
        if min(vnow) < 2.5
            continue
        end
        
        ptnow = pot(tindx==tt);
        tnow = [min(pot(tindx==tt)) max(pot(tindx==tt))];
            
        % get spikes on this trajectory
        snow = [spx(spt>tnow(1) & spt<tnow(2)) spy(spt>tnow(1) & spt<tnow(2))];
        stnow = spt(spt>tnow(1) & spt<tnow(2));
        spnow = spp(spt>tnow(1) & spt<tnow(2));
        
        % convert trajectory into polar coordinates with centroid as origin
        % angle of points
        ptheta = angle2Points(ones(size(pnow)).*ctroi,pnow); % angle between point, centroid and x-axis  

        % distance between point and centroid, point and edge, ratio of the two
        prho = NaN(size(pnow,1),1);
        for pp = 1:length(pnow)
            d1 = sqrt(sum((pnow(pp,:)-ctroi).^2,2)); % distance from point to centroid
            r1 = createRay(ctroi,pnow(pp,:)); % ray from centroid to point
            i1 = intersectRayPolygon(r1,chull); % where ray intersects field edge
            d2 = sqrt(sum((i1-ctroi).^2,2)); % distance to this point
            prho(pp,1) = abs(d1(1))/abs(d2(1)); % ratio of d1 to d2
        end
        
        % convert path data back to cartesian coordinates
        [X,Y] = pol2cart(ptheta,prho);
        pnow2 = [X(:),Y(:)];
        
        % get line between entry/escape
        tline = createLine(pnow(1,:),pnow(end,:));
        
        % angle of this line
        aline = -lineAngle(tline);
        
        % rotate polar coordinates
        ptheta = ptheta + aline;
        
        % get polar coordinate locations of spikes
        stheta = interp1(ptnow,ptheta,stnow);
        srho = interp1(ptnow,prho,stnow);
        
        % accumulate data
        pmat = [pmat; prho ptheta];
        smat = [smat; srho stheta spnow];        

    end
end
if isempty(smat)
    return
end
% perform circular-linear regression
[X,Y] = pol2cart(smat(:,2),smat(:,1));
rdata = getCIRCLINregress(X(:),smat(:,3));
% [beta,R2,beta_ts,R2_ts] = CircularRegression(double(X(:)),double(smat(:,3))); % Thanks Zugaro!
 
% accumulate data
pdata.n_trajectories = ntraj;
pdata.regression_info = rdata;
pdata.regression_points = single([X(:),smat(:,3)]);


if 0
    figure
    subplot(2,3,1)
    im = imagesc(ratemap);
    daspect([1 1 1]);
    axis xy
    title('Map')    

    subplot(2,3,2)
    im = imagesc(mapnow);
    daspect([1 1 1]);
    axis xy
    title('Field')

    subplot(2,3,3)
    plot(pox,poy,'k');
    hold on
    plot(spx,spy,'rx');
    plot(chull(:,1),chull(:,2),'b')
    daspect([1 1 1]);
    axis xy    

    subplot(2,3,4)
    scatter(X,Y,10,smat(:,3));

    subplot(2,3,6)  
    [beta,R2,beta_ts,R2_ts] = CircularRegression(double(X(:)),double(smat(:,3)))

    adata = [X(:) smat(:,3); X(:) smat(:,3)+2*pi; X(:) smat(:,3)+4*pi; X(:) smat(:,3)+6*pi];
    adata(:,2) = rad2deg(adata(:,2));
    plot(adata(:,1),adata(:,2),'k.','MarkerSize',10)

    beta = rad2deg(beta);
    r = refline(beta(1),beta(2));
    set(r,'Color','r','LineWidth',2);
    r = refline(beta(1),beta(2)-360);
    set(r,'Color','r','LineWidth',2);
    r = refline(beta(1),beta(2)+360);
    set(r,'Color','r','LineWidth',2);
    r = refline(beta(1),beta(2)+720);
    set(r,'Color','r','LineWidth',2); 
    r = refline(beta(1),beta(2)+1080);
    set(r,'Color','r','LineWidth',2);     

    ax = gca;
    ax.YLim = [0 720];     
    
    keyboard
end




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




