function [clus] = get_clu_for_klustest(dataformat,config,tets,data_dirs)
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
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
switch dataformat
    case {'Kwikcut'}
        clus = cell(1,size(tets,1));
       
        for tt = 1:size(tets,1)
            electrode_type = tets(tt,2);
            if electrode_type==0 % tetrode
                cutname = [config.cname '_' num2str(tets(tt,1)) '.ntt']; % the cut file name                
            elseif electrode_type==1 % stereotrode
                cutname = [config.cname '_' num2str(tets(tt,1)) '.nst']; % the cut file name                                
            end
            clu = Nlx2MatSpike(cutname,[0 0 1 0 0],0,1,1); 
            clus{tt} = uint8( clu(:) );  

            nclus = numel(unique(clus{tt}));
            disp(sprintf('\t\t...%d total clusters %s',nclus))
        end
        
    case {'Neuralynx'}
        clus = cell(1,size(tets,1));
        for tt = 1:size(tets,1)
            cutname = [config.cname '_' num2str(tets(tt,1)) '.cut']; % the cut file name
            clus{tt} = uint8( getcut(cutname) );  

            nclus = numel(unique(clus{tt}));
            disp(sprintf('\t\t...%d total clusters %s',nclus))
        end
        
    case {'ElePhy'}
        clus = cell(1,max(tets));
        for tt = 1:length(tets)
            for ff = 1:length(data_dirs)
                cutname = [data_dirs{ff} '\' data_dirs{ff} '.mountainsort\output_T' num2str(tets(tt)) '\phy_MS\spike_clusters.npy'];
                clu = readNPY(cutname);

                % ignore noise clusters
                fname = [data_dirs{ff} '\' data_dirs{ff} '.mountainsort\output_T' num2str(tets(tt)) '\phy_MS\cluster_group.tsv'];;
                fid = fopen(fname);
                A1 = textscan(fid,'%s\t%s',1); % column headers
                A2 = textscan(fid,'%d\t%s'); % {clusters} {groups}
                fclose(fid);     
                clu_good = A2{1,1}(ismember(A2{1,2},{'good'}));
                clu(~ismember(clu,clu_good)) = 0;
                clus{tets(tt)} = [clus{tets(tt)}; uint8( clu(:) )];  
            end
            clu_n = unique(clu);
            nclus = length(clu_n);
            disp(sprintf('\t\t...%d total clusters %s',nclus))
        end

end



































