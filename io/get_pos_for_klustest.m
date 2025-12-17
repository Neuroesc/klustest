function [pos,data_intervals,tstart] = get_pos_for_klustest(formats,data_dirs,snames,new_srate,mapset)
% get_pos_for_klustest load position data from .pos files or equivalent
% Data loading function for klustest, given a data format type, and some
% settings, finds .pos files (Tint format) or .nvt files (Neuralynx format) and
% extracts positions, smoothes them, extrapolates missing data, removes jumps in
% the data, computes head direction, angular head velocity and displacement.
%
% USAGE
%
% [pos,data_intervals,tstart] = get_pos_for_klustest(formats,data_dirs,snames,new_srate,mapset)
%
% INPUT
%
% 'formats' - Structure, contains the format type but also information about the tracking LEDs
%
% 'data_dirs' - Data directories, output from get_tets_for_klustest
%
% 'snames' - Session names we want to load, output from get_tets_for_klustest
%
% 'new_srate' - New sample rate we want for the position data, Hz
%
% 'mapset' - Structure, configuration file from klustest
%
% OUTPUT
%
% 'pos' - Table, 'pox','poy','pot','pov','poh','pod','poa'
%       (x, y, time, speed, head direction, displacement direction, angular head
%       velocity)
%
% 'data_intervals' - Nx2, start and end time of each session
%
% 'tstart' - Start time of the recordings, can be added later to get actual time
%
% NOTES
% 1. Phy data format is experimental
%
% 2. Function is not really intended to be used without klustest
% 
% SEE ALSO kwiktint klustest read_rawpos Nlx2MatVT

% HISTORY
%
% version 1.0.0, Release 24/08/16 Initial release, decided to create specific loading functions for klustest
% version 1.0.1, Release 04/04/18 Commenting and formatting for klustest update
% version 1.0.2, Release 29/11/25 Updates for GitHub release
% version 1.0.3, Release 16/12/25 updated filenames for cross platform flexibility
% version 1.0.4, Release 16/12/25 improved comments, not ideal but detailed enough for now
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2025 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    dataformat = formats.pos;
    total_duration = 0;

    switch lower(dataformat)
        case {'kwikcut','neuralynx','phy','kwiktint'}
        % we want to load position data directly from Neuralynx data files
        % i.e. .nvt files
            data_intervals = NaN(length(data_dirs),2);
            pox = [];
            poy = [];
            pot = [];
            pov = [];
            poh = [];
            pod = [];
            poa = [];
            for ff = 1:length(data_dirs) % for every data directory  
                disp(sprintf('\t\t loading'))
                switch dataformat
                    case {'kwiktint'}
                        [~,b,~] = fileparts(data_dirs{ff});
                        h = get_dacq_headers([b '.pos']);
  
                        [led_pos,potn,~] = read_rawpos([b '.pos'],2);
                        potn = potn + total_duration;
                        total_duration = total_duration + h.duration;

                        % extract the coloured LED points
                        nsamps = numel(potn);    
                        tdata = zeros(nsamps,6); 
                        tdata(:,1:2) = [led_pos(:,1,1) led_pos(:,1,2)];
                        tdata(:,3:4) = [led_pos(:,2,1) led_pos(:,2,2)];

                    case {'kwikcut','neuralynx','Neuralynx'}
                        [~,b,~] = fileparts(data_dirs{ff});
                        [potn, ~, ~, ~, targets, ~, ~] = Nlx2MatVT(fullfile(b,'VT1.nvt'), [1 1 1 1 1 1], 1, 1, [] );  
                        potn = potn / 1e6; % convert microseconds to seconds
                        
                        % extract the coloured LED points
                        nsamps = numel(potn);    
                        tdata = zeros(nsamps,6); 
                        for i = 1:nsamps
                            for j = 1:4 % only need to use the first 4 values           
                                target_value = targets(j,i);            
                                binary_target_string = dec2bin(target_value,32);
                    
                                x = bin2dec(binary_target_string(21:32));
                                y = bin2dec(binary_target_string(5:16));            
                    
                                if bin2dec(binary_target_string(2)) % red
                                    tdata(i,1:2) = [x y];
                                elseif bin2dec(binary_target_string(3)) % green
                                    tdata(i,3:4) = [x y];
                                elseif bin2dec(binary_target_string(4)) % blue
                                    tdata(i,5:6) = [x y];
                                end
                            end
                        end

                    case {'phy'}
                        cnames = dir( fullfile(data_dirs{ff},[data_dirs{ff} '*.videoPositionTracking*']) ); % list all cut files corresponding to the given combined name
                        clist = {cnames.name}.'; % get their names                        
                        dat = readTrodesExtractedDataFile( fullfile(data_dirs{ff},clist{1}) );    
                        potn = dat.fields(find(strcmp({dat.fields.name},'time'))).data;
                        potn = double(potn) .* (1/30000); % convert samples to seconds
                        poxn = double(dat.fields(find(strcmp({dat.fields.name},'xloc'))).data);
                        poyn = double(dat.fields(find(strcmp({dat.fields.name},'yloc'))).data);
                        poxn2 = double(dat.fields(find(strcmp({dat.fields.name},'xloc2'))).data);
                        poyn2 = double(dat.fields(find(strcmp({dat.fields.name},'yloc2'))).data);                        
                        % I think this works (if xloc and xloc2 represent 2 LEDs) but the sample data I have has no second LED
                        % also it would probably be better if this was calculated after smoothing the data and removing jumps
                        [pohn,~] = cart2pol(poxn-poxn2,poyn-poyn2); 
                        pohn = rad2deg(pohn);
                end
                data_intervals(ff,:) = [min(potn) max(potn)]; % start and end of this session
                total_time = diff(data_intervals(ff,:));
                disp(sprintf('\b | %.1fs (%.1f mins)',total_time,total_time/60))

                %% post-process position data
                % resample to nrate Hz (actual samples can be irregularly spaced)
                disp(sprintf('\b | resampling'))        
                pot_pp = min(potn) : (1/new_srate) : max(potn);

                % interpolate color data
                tdata2 = zeros(length(pot_pp),6);
                tdata(tdata==0) = NaN;
                for tt = 1:size(tdata,2)
                    tdata2(:,tt) = interp1(potn(:),tdata(:,tt),pot_pp(:),'linear',NaN);
                end
                potn = pot_pp;

                % convert to cm
                if numel(mapset.ppm)==length(data_dirs)
                    ppm = mapset.ppm(ff);
                else
                    ppm = mapset.ppm(1); 
                end        
                tdata2 = tdata2 ./ ppm .* 100;

                % process data
                cindx = reshape((1:size(tdata,2))',2,[])';
                for tt = 1:size(cindx,1)
                    % remove jumps in the data
                    disp(sprintf('\b | jumps'))                
                    pos = [tdata2(:,cindx(tt,1)) tdata2(:,cindx(tt,2))];
                    p2 = circshift(pos,-1,1); % shift data one time point backward
                    ds = sqrt(sum((pos-p2).^2,2)); % distance between these two, gives distance travelled from previous sample at every sample
                    ds(end) = NaN; % last sample cannot be computed
                    zs = ( ds(:) - mean(ds(:),'omitnan') ) / std(ds(:),'omitnan'); % zscore the distances
                    jumps = logical(movmax(zs>mapset.jumpcut,mapset.jumpwindow,'omitnan','Endpoints','shrink')); % smooth this logical vector to try and catch both sides of every jump
                    tdata2(jumps,cindx(tt,1)) = NaN;
                    tdata2(jumps,cindx(tt,2)) = NaN;
                    disp(sprintf('\b (%.2f%% discarded)', sum(jumps)/numel(potn)))                
    
                    % smooth position data using smoothn algorithm (also interpolates NaNs)
                    disp(sprintf('\b | smothing'))                        
                    [datout, ~, ~] = smoothn({tdata2(:,cindx(tt,1)) tdata2(:,cindx(tt,2))},10,'robust');
                    tdata2(:,cindx(tt,1)) = datout{1};
                    tdata2(:,cindx(tt,2)) = datout{2};
                end

                % head direction
                p1 = tdata2(:,cindx(formats.front_led_color,:)); % front LED data
                p2 = tdata2(:,cindx(formats.back_led_color,:)); % back LED data
                p3 = [mean([p1(:,1) p2(:,1)],2) mean([p1(:,2) p2(:,2)],2)]; % average position of the LEDs
    
                % the animal's azimuth is defined by the line from the back LED to the front LED
                % i.e. the vector which points directly 'front-to-back' along the LED array
                [pohn,~] = cart2pol(p1(:,1)-p2(:,1),p1(:,2)-p2(:,2)); 
                pohn = rad2deg(pohn);

                % displacement direction and running speed
                disp(sprintf('\b | speed'))                        
                p3b = circshift(p3,-1,1); % shift data one time point backward
                [podn,ds] = cart2pol(p3(:,1)-p3b(:,1),p3(:,2)-p3b(:,2)); 
                podn = rad2deg(podn);                                
                ds(end) = NaN; % last sample cannot be computed
                wsize = 3;
                dist_sum = movsum(ds,wsize,'omitnan','Endpoints','shrink'); % moving sum of distance over window of wsize samples
                tsize = wsize .* (1/new_srate); % the length (s) of each window
                povn = dist_sum ./ tsize; % speed in cm per second = distance in cm / time in seconds

                % angular head velocity
                pohn(pohn==0) = NaN;
                pahv = rad2deg(angdiff([deg2rad(pohn(:));NaN])) ./ diff([potn(:);NaN]);
                pahv = medfilt1(pahv(:),3,'omitnan');                
                pahv = fillmissing(pahv,'linear');

                % accumulate position data
                pox = [pox; p3(:,1)];
                poy = [poy; p3(:,2)];
                pot = [pot; potn(:)]; 
                pov = [pov; povn(:)]; 
                poh = [poh; pohn(:)];
                pod = [pod; podn(:)];                
                poa = [poa; pahv(:)];

            end
            tstart = min(pot(:));
            pot = pot - tstart; % reset so session starts at t=0    
            pos_array = [pox(:) poy(:) pot(:) pov(:) poh(:) pod(:) poa(:)];    
            pos = array2table(single(pos_array),'VariableNames',{'pox','poy','pot','pov','poh','pod','poa'});
            data_intervals = data_intervals - tstart; % reset so session starts at t=0 
            
        case {'reconstruction'}
        % we want to load position data that was previously reconstructed

            data_intervals = NaN(length(data_dirs),2);
            rpox = []; 
            rpoy = []; 
            rpoz = []; % initialise variables
            gpox = []; 
            gpoy = []; 
            gpoz = []; % initialise variables
            bpox = []; 
            bpoy = []; 
            bpoz = []; % initialise variables
            pov = []; 
            poh = [];
            poa = [];            
            pot = [];
            pod = [];
            for ff = 1:length(data_dirs) % for every data directory  
                % load neuralynx time values
                disp(sprintf('\t\t loading'))
                [potn, ~, ~, ~, ~, ~, ~] = Nlx2MatVT( fullfile(data_dirs{ff},'VT1.nvt'), [1 1 1 1 1 1], 1, 1, [] );  
                potn = potn ./ 1e6; % convert microseconds to seconds
                data_intervals(ff,:) = [min(potn) max(potn)]; % start and end of this session
                total_time = diff(data_intervals(ff,:));
                disp(sprintf('\b | %.1fs (%.1f mins)',total_time,total_time/60))   
                
                % load 3D trajectory data
                fcheck = fullfile(pwd,'3Dreconstruction',[snames{ff} '_merged.txt']);
                if exist(fcheck,'file') % check it exists
                    cdata = readtable(fcheck,'Delimiter','\t');
                else
                    error('Predetermined 3D trajectory not found: %s\n ...exiting',fcheck)
                end
               
                % resample to nrate Hz
                disp(sprintf('\b | resampling'))        
                potn = double( cdata.pot(:) ./ 1e6 );
                pot_pp = min(potn) : (1/new_srate) : max(potn);    
                rpoxn = interp1(potn,double( cdata.rpox(:) ),pot_pp,'linear',NaN);
                rpoyn = interp1(potn,double( -cdata.rpoy(:) ),pot_pp,'linear',NaN);
                rpozn = interp1(potn,double( -cdata.rpoz(:) ),pot_pp,'linear',NaN);
                rpox = [rpox(:); rpoxn(:)];
                rpoy = [rpoy(:); rpoyn(:)];
                rpoz = [rpoz(:); rpozn(:)];
                
                gpoxn = interp1(potn,double( cdata.gpox(:) ),pot_pp,'linear',NaN);
                gpoyn = interp1(potn,double( -cdata.gpoy(:) ),pot_pp,'linear',NaN);
                gpozn = interp1(potn,double( -cdata.gpoz(:) ),pot_pp,'linear',NaN);                
                gpox = [gpox(:); gpoxn(:)];
                gpoy = [gpoy(:); gpoyn(:)];
                gpoz = [gpoz(:); gpozn(:)]; 
                
                bpoxn = interp1(potn,double( cdata.bpox(:) ),pot_pp,'linear',NaN);
                bpoyn = interp1(potn,double( -cdata.bpoy(:) ),pot_pp,'linear',NaN);
                bpozn = interp1(potn,double( -cdata.bpoz(:) ),pot_pp,'linear',NaN); 
                bpox = [bpox(:); bpoxn(:)];
                bpoy = [bpoy(:); bpoyn(:)];
                bpoz = [bpoz(:); bpozn(:)];                
                potn = pot_pp; 
                pot = [pot(:); potn(:)];

                % calculate roll pitch yaw
                disp(sprintf('\b | 3D HD'))                                                        
                P0 = [rpoxn(:) rpoyn(:) rpozn(:)]; % red position
                P1 = [bpoxn(:) bpoyn(:) bpozn(:)]; % blue position
                P2 = [gpoxn(:) gpoyn(:) gpozn(:)]; % green position  
                mpos = [mean([rpoxn(:) gpoxn(:) bpoxn(:)],2,'omitnan') mean([rpoyn(:) gpoyn(:) bpoyn(:)],2,'omitnan') mean([rpozn(:) gpozn(:) bpozn(:)],2,'omitnan')]; % mean position
                % the red, green and blue LEDs form a plane, this next line will find the surface normal of this
                % i.e. the vector which points directly 'up' perpendicular from the LED array
                % we can use this to correct the position data and get the animal's head position
                % note that if the surface normal z coordinate is below the medium point the animal must be inverted
                led_normal = normalize( cross(P0-P1, P0-P2),2,'norm' ); 
                
                % correct the mean position point so it is at the animal's head
                hpos = mpos - (led_normal*mapset.drive_height_mm); % hpos = head position
                
                % the red, mean position and surface normal form a plane, this next line will find the surface normal of this
                % i.e. the vector which points directly 'left-to-right' along the LED array      
                % we can use this to calculate the animal's head roll
                roll_normal = normalize( cross(P0-mpos, P0-(mpos+led_normal)),2,'norm' );  
                % negative roll means left ear is lower than the right
                
                % the animal's azimuth is defined by the line from the mean position to the red LED
                % i.e. the vector which points directly 'front-to-back' along the LED array
                % we can use this to calculate the animal's yaw and pitch
                nose_normal = normalize( P0-mpos,2,'norm' );  
                % yaw: [-y 0] = -90, [+y 0] = 90, [+x 0] = 0, [-x 0] = 180
                
                % yaw, pitch and roll from these
                [yaw,pitch,~] = cart2sph(nose_normal(:,1),nose_normal(:,2),nose_normal(:,3));
                [~,roll,~] = cart2sph(roll_normal(:,1),roll_normal(:,2),roll_normal(:,3));
                roll(led_normal(:,3)<0) = roll(led_normal(:,3)<0)+pi; % if the animal is inverted add 180 degrees to roll
                head_angles = array2table([hpos yaw pitch roll],'VariableNames',{'meanRx','meanRy','meanRz','yaw','pitch','roll'});                
                poh = [poh; head_angles];                
                
                % running speed
                disp(sprintf('\b | speed'))                        
                pos = [head_angles.meanRx(:) head_angles.meanRy(:) head_angles.meanRz(:)];
                p2 = circshift(pos,-1,1); % shift data one time point backward
                [podn,~] = cart2pol(pos(:,1)-p2(:,1),pos(:,2)-p2(:,2)); 
                ds = sqrt(sum((pos-p2).^2,2)); % distance between these two, gives distance travelled from previous sample at every sample
                ds(end) = NaN; % last sample cannot be computed
                podn(end) = NaN;
                pod = [pod; podn(:)];
                wsize = 3;
                dist_sum = movsum(ds,wsize,'omitnan','Endpoints','shrink'); % moving sum of distance over window of wsize samples
                tsize = wsize .* (1/new_srate); % the length (s) of each window
                povn = dist_sum ./ tsize; % speed in cm per second = distance in cm / time in seconds     
                pov = [pov(:); povn(:)];                
                
                % angular head velocity
                pahv = diff(rad2deg(head_angles.yaw(:))) ./ diff(potn(:));
                pahv = medfilt1(pahv(:),3,'omitnan');
                poa = [poa; pahv(1); pahv(:)];
                
            end
            tstart = min(pot(:));
            pot = pot(:) - tstart; % reset so session starts at t=0 
            data_intervals = data_intervals - tstart; % reset so session starts at t=0   
            ppm = 1000;
            pos_array = [poh.meanRx(:)./ppm.*100 poh.meanRy(:)./ppm.*100 poh.meanRz(:)./ppm.*100 pot(:) pov(:) poa(:) rad2deg(poh.yaw(:)) rad2deg(pod(:)) rad2deg(poh.pitch(:)) rad2deg(poh.roll(:))];    
            pos = array2table(single(pos_array),'VariableNames',{'pox','poy','poz','pot','pov','poa','poh','pod','pitch','roll'});                
  
    end






































