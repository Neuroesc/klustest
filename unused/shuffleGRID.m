function [jtimes,tshuff,sos] = shuffleSPIKES(duration,itimes,mint,pos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function randomly shuffles spike trains by a random time value (+ve or -ve)
%   If the position data (x,y,t) for the session are also given 
%   jtimes = shuffleSPIKES(duration,itimes,mint)
%
%%%%%%%% Inputs
%   duration = (s) the maximum length of the recording session
%   itimes = (s) the spike times
%   mint = (default = 20s) the minimum time we want to shuffle the train by
%   pos (optional) = the position data [x y t] coordinates
%
%%%%%%%% Outputs
%   jtimes = the new, shuffled spike train
%   tshuff = the time the train was shuffled by
%   sos (if pos is given) = the new spike [x y t] coordinates that correspond to the observed data, calculated using a nearest neighbour approach
%
%%%%%%%% Comments
%   10/08/16 created 
%   11/08/16 actually got it working
%   © Roddy Grieves: rmgrieves@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform shuffle
                %% Shuffled grid analysis
                disp('shuffling')
                reps = 500;
                g_shuff = NaN(reps,1);
                parfor i = 1:reps
                    [jtimes,~,posout] = shuffleSPIKES(duration,spt,20,[pox poy pot]);
                    spt_new = posout(:,3);
                    spx_new = posout(:,1);
                    spy_new = posout(:,2);
                    
                    % make time/dwell/rate maps
                    [dwellmap,~] = mapDATA(pox,poy,map_limits,bin_size,pixel_ratio);
                    dwellmap = dwellmap .* pos_tb;
                    dwellmap(dwellmap < min_dwell) = 0;
                    dwellmap = imgaussfilt(dwellmap,map_sigma);
                    [spikemap,~] = mapDATA(spx_new,spy_new,map_limits,bin_size,pixel_ratio);
                    spikemap = imgaussfilt(spikemap,map_sigma);
                    ratemap = spikemap ./ dwellmap;
                    ratemap(dwellmap == 0) = NaN;

                    % make and analyse autocorrelation
                    automap = GridAutoCorr(ratemap);
                    [grid_score,~,~,~,~] = GridAnalysis(automap,bin_size);
                    g_shuff(i,1) = grid_score;
                end % for i = 1:reps
                sdata.(tet_s).(clu_s).grid_score_shuffled_distribution = g_shuff; % add data to structure
                disp('done')

