function pos = kamkwik(pos,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
%FUNCTION  short desc.
%
% USAGE:
%           out = template(in) process with default settings
%
%           out = template(in,optional1) process using optional argument 1
%
%           out = template(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%           of the process
%
%           Parameters include:
%
%           'param1'          -   Scalar value, parameter to do something
%
%           'param2'          -   Scalar value, parameter to do something
%
% EXAMPLES:
%
%           % run function using default values
%           out = template(in,varargin)
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2019 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INPUT ARGUMENTS CHECK
%% Prepare default settings
    def_cshots              = 1;        
    
%% Parse inputs
    p = inputParser;
    addRequired(p,'pos',@(x) ~isempty(x) && isnumeric(x));  
    addParameter(p,'cshots',def_cshots,@(x) isnumeric(x) && isscalar(x));   
    parse(p,in,varargin{:});

%% Retrieve parameters 
    config = p.Results;
    
%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
% get the 




















































































