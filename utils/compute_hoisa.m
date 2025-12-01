function [f,xi,e,edg] = compute_hoisa(spt1,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS
% FUNCTION  short desc.
% long description
%
% log-sampled ISIHs: https://pmc.ncbi.nlm.nih.gov/articles/PMC6772912/
% HOISA: https://link.springer.com/chapter/10.1007/978-3-642-20853-9_32
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
%       % random spike trains
%       s1 = cumsum(rand([n1,1])) .* 1000; % ms
%       s2 = cumsum(rand([n2,1])) .* 1000; % ms
%
%       % poisson spike train, 100ms burst and a monosynaptic cell, +20ms offset
%       s1 = cumsum(poissrnd(100,[n1,1])); % ms
%       s2 = s1+abs(normrnd(20,10,[n1,1]));
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 19/12/24 Initial release
%
% AUTHOR 
% Roddy Grieves
% University of Glasgow, Sir James Black Building
% Neuroethology and Spatial Cognition Lab
% eMail: roddy.grieves@glasgow.ac.uk
% Copyright 2024 Roddy Grieves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION BODY
    % Parse inputs
    p = inputParser;
    addRequired(p,'spt1',@(x) isnumeric(x));  
    addParameter(p,'spt2',[],@(x) isnumeric(x));   
    addParameter(p,'binsize',1,@(x) isnumeric(x) && isscalar(x));   
    addParameter(p,'window',500,@(x) isnumeric(x) && isscalar(x)); 
    addParameter(p,'max_diffs',1e8,@(x) isnumeric(x));   
    addParameter(p,'min_diffs',1e4,@(x) isnumeric(x));       
    addParameter(p,'method','hoisa',@(x) any(strcmp(x,{'hoisa','isih'})));       
    parse(p,spt1,varargin{:});
    config = p.Results;

    % prepare data
    s1 = config.spt1 .* 1000; % convert to ms
    s2 = config.spt2 .* 1000; % convert to ms
    % s1 = cumsum(poissrnd(100,[5000,1]));
    % s1 = cumsum(poissrnd(100,[5000,1])+normrnd(0,15,[5000,1]));

    % prepare bins
    edg = 0 : config.binsize : config.window;
    edg = unique([-edg edg]);
    xi = movmean(edg,2,'Endpoints','discard');

    n1 = numel(s1);
    n2 = numel(s2);
    tot_n = n1*n2;
    
    % generate outputs
    e = false;
    auto = isequal(s1,s2);
    switch config.method
        case {'hoisa','isih'} % Higher Order Interspike Autocorrelation or Interspike interval histogram
            if tot_n<config.max_diffs
                tdiff = s2 - s1'; % logical expansion of difference
                idx = tdiff>-config.window & tdiff<config.window;
                if auto
                    idx(1 : n2+1 : n2^2) = false;
                end  
                if strcmp(config.method,'isih')
                    tdiff(~idx | tdiff<0) = NaN;
                    m = min(tdiff,[],2,'omitmissing');
                    f = histcounts(m,edg,'Normalization','count');
                else
                    f = histcounts(tdiff(idx),edg,'Normalization','count');                
                end
                f = f ./ sum(f(:),'all','omitmissing');

            else
                maxcorr = min([config.min_diffs n1]); % limit the number of comparisons to min_comparisons
                e = maxcorr==config.min_diffs;
                rindx = randperm(n1,maxcorr);
                fs = NaN(length(rindx),length(xi));
                for ff = 1:length(rindx)
                    rej_spike = false(size(s2));
                    if auto
                        rej_spike(rindx(ff)) = true;
                    end
                    tdiff = s2(s2>s1(rindx(ff))-config.window & s2<s1(rindx(ff))+config.window & ~rej_spike) - s1(rindx(ff));
                    if strcmp(config.method,'isih')
                        tdiff(tdiff<0) = NaN;
                        m = min(tdiff,[],'all','omitmissing');
                        fs(ff,:) = histcounts(m,edg,'Normalization','count');
                    else
                        fs(ff,:) = histcounts(tdiff(:),edg,'Normalization','count');
                    end
                end
                f = sum(fs,1,'omitmissing');
                f = f ./ sum(f(:),'all','omitmissing');
            end

            f(isnan(f)) = 0;            
            if (auto || strcmp(config.method,'isih'))
                f = [fliplr(f(xi>=0)) f(xi>=0)];
                xi = unique([fliplr(-xi) xi]);
            end
    end

























