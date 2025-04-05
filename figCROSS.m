function figCROSS(sdata,varargin)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% figCROSS A figure showing the cross- spike correlogram for all clusters on a tetrode
% This function cross-correlates spike trains for every cluster pair on a tetrode
% and plots them so that overlapping clusters can be identified. This can also be
% used to determine if spikes from one cluster consistently occur before/after another.
%
% In the cross-correlograms, a peak after zero would mean spikes in the first spike train occur 
% frequently right after the spikes in the second train. In this figure train 1 = cell # shown at 
% the top, train 2 = cell # shown at the side
%
% The waveforms shown are the mean waveform for the cluster taken from the channel with the
% highest amplitude.
%
% USAGE:
%       figCROSS(sdata) process with default settings
% 
%       figCROSS(sdata___,Name,Value,...) process with Name-Value pairs used to control aspects 
%       of the process
% 
%       Parameters include:
% 
%       'tet'          -   (default = 1) Scalar, which tetrode to analyse
% 
%       'figfile'      -   (default = pwd) String, directory, where to save the figure
%                          This doesn't include the figure name, just the path to the directory where it should be saved
% 
%       'fig_vis'      -   (default = 'off') String, 'on' = show figure when plotting, 'off' = hide figure
% 
%       'corr_window'  -   (default = 60) Scalar, time in ms over which correlations should be plotted (i.e. min/max on x-axis)
% 
%       'corr_binsize' -   (default = 0.75) Scalar, duration in ms of bins to use (i.e. width of bars)
%
% INPUT:
%       sdata    - sdata table provided by klustest or loaded from a file
%
% EXAMPLES:
%       % run function using default values
%       figCROSS(sdata,'tet',1,'figfile',pwd,'fig_vis',0)
%
% See also: klustest spikeINTERVALS_v2

% HISTORY:
% version 1.0.0, Release 31/03/17 Initial release
% version 2.0.0, Release 20/01/22 Updated to work with table format sdata
% version 2.0.1, Release 20/01/22 Added reflection of symmetrical cross-correlations to save time
% version 2.1.0, Release 20/01/22 Added waveforms to help identify overlapping clusters
% version 2.1.1, Release 20/01/22 Updated comments and varagin input
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
    addRequired(p,'sdata',@(x) ~isempty(x) );  
    addParameter(p,'tet',1,@(x) isnumeric(x) && isscalar(x));  
    addParameter(p,'figfile',pwd,@(x) ischar(x));  
    addParameter(p,'fig_vis','off',@(x) ischar(x));  
    addParameter(p,'corr_window',60,@(x) isnumeric(x) && isscalar(x));  
    addParameter(p,'corr_binsize',0.75,@(x) isnumeric(x) && isscalar(x));  
    
    parse(p,sdata,varargin{:});
    config = p.Results;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% Initial variables
    sdata = config.sdata;
    clusters = sdata.clu(sdata.tet==config.tet); % list of available clusters
    clus_count = numel(clusters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot figure
    fig_corr = figure('visible',config.fig_vis,'Position',[100, 100, 1200, 850]);
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    fig_hor = clus_count; % how many plots wide should it be
    fig_ver = clus_count; % how many plots tall should it be
    fspac = 0.005; % the spacing around the plots, on all sides
    fpadd = 0.005; % the spacing around the plots, on all sides, this takes more space than fspac though
    fmarg = 0.05; % the margins around the plots, at the edge of the figure

%% >>>>>>>>>>>>>>>>>>>> plot observed cluster auto-correlograms diagonally
    cols = inferno(clus_count);
    for cc = 1:clus_count
        idx = find( sdata.tet==config.tet & sdata.clu==clusters(cc) );
        cspt = sdata.spike_times{idx};

        sindx = sub2ind([clus_count,clus_count],cc,cc);    
        subaxis(fig_ver,fig_hor,sindx,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
            [tisi,sisi] = spikeINTERVALS_v2(cspt,config.corr_window,config.corr_binsize);  
            bar(tisi,sisi,1,'k')
            ax = gca;
            ax.XLim = [-config.corr_window,config.corr_window];  
            ax.FontSize = 8;    
            if cc==clus_count
                xlabel('Time lag (ms)')
            else
                ax.XTick = []; % change Xtick locations to these values            
            end       
            if cc==1
                ylabel('Spikes')
                text(0.5,1.2,sprintf('Cell %d',cc),'HorizontalAl','center','VerticalAl','middle','Units','normalized','FontSize',14);
                text(-0.2,0.5,sprintf('Cell %d',cc),'HorizontalAl','center','VerticalAl','middle','Units','normalized','rotation',90,'FontSize',14);            
            else
                ax.YColor = 'none';
            end
            box off

            % plot mini-waveform in a sub-plot
            ax.Units = 'pixels';
            ax2pos = [ax.Position(1)+(ax.Position(3)*0.02), ax.Position(2)+ax.Position(4)*0.66, ax.Position(3)*0.33, ax.Position(4)*0.33];
            ax2 = axes(fig_corr,'Units','pixels','Position',ax2pos);
                [~,max_wav] = max(sdata.waveform_max(idx,:));
                wav = sdata.waveform_mean{idx,max_wav};
                wavtime = -200:20:780;
                plot(wavtime,wav,'Color',cols(cc,:));
                axis off    
    end

%% >>>>>>>>>>>>>>>>>>>> plot cluster cross-correlations
    pairs = nchoosek(1:clus_count,2); % every possible pair of clusters, excluding auto-correlations
    for pp = 1:length(pairs(:,1))
        pnow = pairs(pp,:);
        idx1 = find( sdata.tet==config.tet & sdata.clu==clusters(pnow(1)),1 );
        idx2 = find( sdata.tet==config.tet & sdata.clu==clusters(pnow(2)),1 );    
        spt1 = sdata.spike_times{ idx1 };
        spt2 = sdata.spike_times{ idx2 };
        [tisi,sisi] = spikeINTERVALS_v2(spt1,config.corr_window,config.corr_binsize,spt2);  

        sindx = sub2ind([clus_count,clus_count],pnow(1),pnow(2));
        subaxis(fig_ver,fig_hor,sindx,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
            bar(tisi,sisi,1,'b');
            ax = gca;
            ax.XLim = [-config.corr_window,config.corr_window];        
            if pnow(2)==clus_count
                xlabel('Time lag (ms)')
            else
                ax.XTick = []; % change Xtick locations to these values            
            end       
            if pnow(1)==1
                ylabel('Spikes')
                text(-0.2,0.5,sprintf('Cell %d',clusters(pnow(2))),'HorizontalAl','center','VerticalAl','middle','Units','normalized','rotation',90,'FontSize',14);            
            else
                ax.YColor = 'none';
            end        
            ax.FontSize = 8;
            box off

            % plot mini-waveform
            ax.Units = 'pixels';
            ax2pos = [ax.Position(1)+(ax.Position(3)*0.02), ax.Position(2)+ax.Position(4)*0.66, ax.Position(3)*0.33, ax.Position(4)*0.33];
            ax2 = axes(fig_corr,'Units','pixels','Position',ax2pos);
                [~,max_wav1] = max(sdata.waveform_max(idx1,:));
                wav1 = sdata.waveform_mean{idx1,max_wav1};
                [~,max_wav2] = max(sdata.waveform_max(idx2,:));
                wav2 = sdata.waveform_mean{idx2,max_wav2};
                wavtime = -200:20:780;
                plot(wavtime,wav1,'Color',cols(pnow(1),:)); hold on;
                plot(wavtime,wav2,'Color',cols(pnow(2),:));           
                axis off        

        sindx = sub2ind([clus_count,clus_count],pnow(2),pnow(1));
        subaxis(fig_ver,fig_hor,sindx,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
            % to save time the cross-correlation for the opposite direction (i.e. cell 1 vs cell 2 instead of cell 2 vs cell 1)
            % is simply plotted here reflected around the y-axis
            bar(-tisi,sisi,1,'b');
            ax = gca;
            ax.XLim = [-config.corr_window,config.corr_window];        
            ax.YTick = []; % change Xtick locations to these values
            ax.YColor = 'none';
            ax.FontSize = 8;
            ax.XTick = []; % change Xtick locations to these values            
            box off
            if pnow(1)==1
                text(0.5,1.2,sprintf('Cell %d',clusters(pnow(2))),'HorizontalAl','center','VerticalAl','middle','Units','normalized','FontSize',14);
            end  

            % plot mini-waveform in a sub-plot
            ax.Units = 'pixels';
            ax2pos = [ax.Position(1)+(ax.Position(3)*0.02), ax.Position(2)+ax.Position(4)*0.66, ax.Position(3)*0.33, ax.Position(4)*0.33];
            ax2 = axes(fig_corr,'Units','pixels','Position',ax2pos);
                plot(wavtime,wav1,'Color',cols(pnow(1),:)); hold on;
                plot(wavtime,wav2,'Color',cols(pnow(2),:));           
                axis off       
    end

    % saveas is faster than print
    saveas(fig_corr,[config.figfile 'electrode_' num2str(config.tet) '_crosscorr.png'],'png')
    close(fig_corr);     

































