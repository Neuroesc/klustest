function alog = analysis_log(analysis,status,varargin)
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
    addOptional(p,'analysis',cell(1,1),@(x) iscell(x)); 
    addOptional(p,'status',1,@(x) isnumeric(x)); 
    addParameter(p,'notes',{'na'},@(x) iscell(x));     
    addParameter(p,'fname',[pwd '\_analysis_log.txt'],@(x) isstring(x) || ischar(x));  
    addParameter(p,'version',{'v0.0'},@(x) iscell(x));      
    parse(p,analysis,status,varargin{:});
    config = p.Results;
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    if exist(config.fname,'file')
        alog = readtable(config.fname,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
        if isempty(config.analysis)
            return
        end
        
        alog_now = table;
        alog_now.analysis = config.analysis;
        alog_now.status = config.status;
        alog_now.date = datestr(now,30);        
        alog_now.notes = config.notes;
        alog_now.version = config.version;
        
        idx = find(strcmp(alog.analysis,config.analysis));
        if isempty(idx)
            alog_now.n = 1;                
            alog = [alog; alog_now];
        else
            alog_now.n = alog.n(idx)+1;                           
            alog(idx,:) = alog_now(1,:);
        end

        writetable(alog,config.fname,'FileType','text','WriteVariableNames',1,'WriteMode','overwrite','Delimiter','\t');        
    else
        alog = table;
        alog.analysis = config.analysis;
        alog.status = config.status;
        alog.date = datestr(now,30);
        alog.notes = config.notes;
        alog.version = config.version;   
        alog.n = 1;   
        
        writetable(alog,config.fname,'FileType','text','WriteVariableNames',1,'WriteMode','overwrite','Delimiter','\t');
    end






























