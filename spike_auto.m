function [acg,xi,estimated] = spike_auto(spt1,varargin)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  short desc.
% long description
%
% USAGE:
%       [out] = template(in) process with default settings
% 
%       [out] = template(in,optional1) process using optional argument 1
% 
%       [out] = template(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%       of the process
% 
%       Parameters include:
% 
%       'param1'          -   (default = X) Scalar value, parameter to do something
% 
%       'param2'          -   (default = X) Scalar value, parameter to do something
% 
% INPUT:
%       in    - input as a vector
%
% OUTPUT:
%       out   - output as a vector
%
% EXAMPLES:
%       % run function using default values
%       out = template(in,varargin)
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    p = inputParser;
    addRequired(p,'spt1',@(x) isnumeric(x));  
    addParameter(p,'spt2',[],@(x) isnumeric(x));   
    addParameter(p,'bin_size',1,@(x) isnumeric(x) && isscalar(x));   
    addParameter(p,'win_size',500,@(x) isnumeric(x) && isscalar(x)); 
    addParameter(p,'log',0,@(x) isnumeric(x) && isscalar(x)); 
    addParameter(p,'method','correlogram',@(x) any(strcmp(x,{'correlogram','isi','j_isi'})));       
    parse(p,spt1,varargin{:});
    config = p.Results;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
estimated = false;
switch config.method
%% >>>>>>>>>>>>>>>>>>>> Correlogram
    case {'correlogram'}
        bsiz = config.bin_size; % ms
        window = config.win_size; % ms
        edg = 0 : bsiz : window; % histogram bin edges
    
        mode = 'cross';
        if isempty(config.spt2) % if this is an autocorrelation
            config.spt2 = config.spt1;
            mode = 'auto';
        end
        
        if isempty(config.spt1) % if there are no spikes
            xi = movmean(edg,2,'EndPoints','discard'); % bin centers for plotting
            acg = zeros(size(xi)); % return an empty histogram 
            return
        end    
        
        spt1 = double( config.spt1(:) ) .* 1000; % convert to ms
        spt2 = double( config.spt2(:) ) .* 1000; % convert to ms
       
        if (numel(spt1)*numel(spt2))<2e6 % if there are not too many spikes
            tdiff = spt1 - spt2'; % logical expansion of difference
            if strcmp(mode,'auto')
                % if this is an autocorrelation, remove differences between spikes
                % and themselves
                idx = 1 : length(spt1)+1 : length(spt1)^2;
                tdiff(idx) = NaN; % remove the current spike
            end  
            f = histcounts(tdiff(:),edg); % bin spike distances

        else % if there are many spikes compute differences one element at a time
            estimated = true;

            f = zeros(1,length(edg)-1);
            maxcorr = min([1e3 numel(spt1)]); % limit the number of comparisons to 1e3
            rng(0);
            rindx = randperm(numel(spt1),maxcorr);
            sidx = true(size(spt2));                    
            for ii = 1:maxcorr
                if strcmp(mode,'auto') % remove the current spike
                    sidx(rindx(ii)) = false;
                end
                tindx = spt2>(spt1(rindx(ii))-window-1) & spt2<(spt1(rindx(ii))+window+1);
                tdiff = (ones(sum(tindx & sidx),1).*spt1(rindx(ii))) - spt2(tindx & sidx); % differences
                fn = histcounts(abs(tdiff(:)),edg); % bin spike distances
                f = f + fn;
                if strcmp(mode,'auto')
                    sidx(rindx(ii)) = true; % restore for next loop
                end               
            end
        end
    
        % reflect values around zero
        f = [flipud(f(:)); f(:)];
        edg = movmean(edg,2,'EndPoints','discard');
        xi = [flipud(-edg(:)); edg(:)];
        acg = f ./ sum(f(:)) ./ bsiz; % normalize to unit interal (pdf)
    
%% >>>>>>>>>>>>>>>>>>>> ISI     
    case {'isi'} % interspike interval
        bsiz = config.bin_size; % ms
        window = config.win_size; % ms
        edg = 0 : bsiz : window; % histogram bin edges
        xi = movmean(edg,2,'EndPoints','discard');

        if isempty(config.spt1) % if there are no spikes
            acg = zeros(size(xi)); % return an empty histogram 
            return
        end  

        mode = 'cross';
        if isempty(config.spt2) % if this is an autocorrelation
            config.spt2 = config.spt1;
            mode = 'auto';
        end  
        
        spt1 = double( config.spt1(:) ) .* 1000; % convert to ms
        spt2 = double( config.spt2(:) ) .* 1000; % convert to ms

        if strcmp(mode,'auto')
            tdiff = diff(spt1); % inter-spike intervals
            if config.log
                edg = log(edg);
                tdiff = log(tdiff);
            end
            f = histcounts(tdiff(:),edg,'Normalization','pdf'); % bin spike distances

        else
            if (numel(spt1)*numel(spt2))<2e6 % if there are not too many spikes
                tdiff = spt1 - spt2'; % logical expansion of difference
    
                tdiff1 = tdiff;
                tdiff1(tdiff>0) = NaN; % keep only future distances
                cross_diff = min(abs(tdiff1),[],2,'omitnan'); % distance from spikes in spt1 to closest next spike in spt2
    
                f = histcounts(cross_diff(:),edg,'Normalization','pdf'); % bin spike distances
    
            else % if there are many spikes compute differences one element at a time
                estimated = true;
    
                maxcorr = min([1e3 numel(spt1)]); % limit the number of comparisons to 1e3 or the number of spikes in spt1, whichever is smallest
                rng(0); % make the result reproducible
                rindx = randperm(numel(spt1),maxcorr); % sampling without replacement
                vals = NaN(maxcorr,1);
                for ii = 1:maxcorr
                    tindx = spt2>=spt1(rindx(ii)) & spt2<spt1(rindx(ii))+window+1;
                    tdiff = spt1(rindx(ii)) - spt2(tindx); % differences
                    if ~isempty(tdiff)
                        vals(ii) = min(abs(tdiff(:)));
                    end
                end
                f = histcounts(vals(:),edg,'Normalization','pdf'); % bin spike distances
                
            end
        end

        acg = [flipud(f(:)); f(:)];
        xi = [flipud(-xi(:)); xi(:)];

end





























