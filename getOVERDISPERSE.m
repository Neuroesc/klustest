function [z,overdispersion,r] = getOVERDISPERSE(mpos,ppot,pspt,ratemap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getOVERDISPERSE  calculate overdispersion as in Fenton et al. (2010)
% Takes position data scaled to match firing rate map, time of position data, time of spikes and firing rate map
% to calculate the overdispersion
% In this case z measures the deviation of observed discharge from expected in standard deviation units, overdispersion 
% in turn is the variance of the z distribution for a set of passes. 
%
% USAGE:
%         [z,overdispersion,r] = getOVERDISPERSE(mpos,pot,spt,ratemap)
%
% INPUT:
%         mpos - position data in map coordinates [x,y]
%         pot - time corresponding to position data
%         spt - spike times
%         ratemap - firing rate map of the session
%
% OUTPUT:
%    z - deviation of observed discharge (firing rate in 5s bins) from expected (firing rate in map) in standard deviation units
%    overdispersion - the variance of the z distribution
%    r - correlation between observed firing rate and expected
%
% See also: klustest accumarray

% HISTORY:
% version 1.0.0, Release 08/04/19 Initial release, contained here instead of sitting in klustest
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with input variables
    narginchk(4,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION BODY
    % Taken from: Fenton et al. (2010) Attention-like modulation of hippocampus place cell discharge
    OD_bin_size = 5; % time window of overdispersion analysis

    % In the first method, the entire session was divided into 5-sec intervals. For each interval we calculated the expected number of spikes, exp, as:
    % [equation in paper]
    % where ri is the time-averaged rate at location i, and ti is the time spent in location i during the pass. Only intervals during which exp ? 5.0 AP were 
    % used to calculate overdispersion since the overall firing rate of place cells is ~1.0 AP/sec.
    % convert position coordinates to map coordinates and extract firing rate values from map
    poxmap = mpos(:,1);
    poymap = mpos(:,2);
    exp_frate = ratemap(sub2ind(size(ratemap),round(poymap),round(poxmap)));

    % bin spikes into OD_bin_size second windows
    edg = min(ppot):OD_bin_size:max(ppot);
    [~,~,bindex] = histcounts(ppot,edg);

    % bin observed spikes into the same windows
    [obs_spikes,~,~] = histcounts(pspt,edg);
    bindex(bindex==0) = nanmax(bindex(:))+1; % some spikes fall outside our maximum bin, we just ignore these

    % calculate the expected number of spikes for each time window, based on the firing rate map
    exp_spikes = accumarray(bindex,exp_frate.*(1/50))';     
    exp_spikes(end) = []; % remove the maximum bin we made before

    % For each selected 5-sec interval, we then calculated z, the normalized standard deviation of obs, the observed number of spikes as:
    % [equation in paper]
    % calculate z, the overdispersion metric
    z = (obs_spikes - exp_spikes) ./ sqrt(exp_spikes);
    z(exp_spikes<5) = [];
    overdispersion = nanstd(z)^2;
    r = corr(obs_spikes(:),exp_spikes(:),'type','Spearman','rows','pairwise');
    % z measures the deviation of observed discharge from expected in standard deviation units. Overdispersion in turn is the variance of 
    % the z distribution for a set of passes. The outcome of this calculation was found to be indistinguishable from the somewhat different 
    % method previously used (Fenton and Muller, 1998).
    % z should have values somewhere between 2 and 5, there is no real upper bound though

    % allocate NaN values to z and overdispersion instead of leaving them empty
    if isempty(z)
        z = NaN;
    end
    if isempty(overdispersion)
        overdispersion = NaN;
    end                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        









