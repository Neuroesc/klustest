function [sdata,pdata] = addCELLTYPE(varargin)
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
    def_sdata              = [];        
    def_pdata              = [];        
    def_cname              = 'kwiktint';        
    def_oldct              = 1;        
    
%% Parse inputs
    p = inputParser;
    addParameter(p,'sdata',def_sdata,@(x) ~isempty(x) && istable(x)); 
    addParameter(p,'pdata',def_pdata,@(x) ~isempty(x) && isstruct(x)); 
    addParameter(p,'cname',def_cname,@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    addParameter(p,'oldct',def_oldct,@(x) isnumeric(x) && isscalar(x));   
    parse(p,varargin{:});

%% Retrieve parameters 
    config = p.Results;
    
%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% Get sdata and pdata
    sdname = [pwd '\klustest\' config.cname '\' config.cname 'sdata.mat']; 
    pdname = [pwd '\klustest\' config.cname '\' config.cname 'pdata.mat']; 
    if isempty(config.sdata) | isempty(config.pdata)
        load(sdname,'sdata');
        load(pdname,'pdata');
    else
        sdata = config.sdata;
        pdata = config.pdata;
    end

    % get a list of all the clusters and parts available
    cname = pdata.combined_name;
    cnames = unique(sdata.uci(:));
    snames = unique(sdata.part(:));
    nsess = length(snames);

%% if we just want to combine previous cell type data with a new sdata (i.e. load previous cell IDs and use them)
    ctname = [pwd '\klustest\' config.cname '\' config.cname '_cell_types.mat']; 
    if exist(ctname,'file') && config.oldct
        load(ctname,'faux_sdata');
        if all(ismember({'auto_cid_name','auto_cid_type','auto_cid_info','manual_cid_name','manual_cid'},sdata.Properties.VariableNames))
            return
        end

        sdata = addToTable(sdata,2,'auto_cid_info',NaN(1,5),'auto_cid_type',NaN(1,1),'auto_cid_name',cell(1,1),'manual_cid',NaN(1,1),'manual_cid_name',cell(1,1)); % [table in, method 1=add 2=preallocate, variable name, variable value]    
        sdata.auto_cid_info(:) = faux_sdata.auto_cid_info(:);
        sdata.auto_cid_type(:) = faux_sdata.auto_cid_type(:);
        sdata.auto_cid_name(:) = faux_sdata.auto_cid_name(:);
        sdata.manual_cid(:) = faux_sdata.manual_cid(:);
        sdata.manual_cid_name(:) = faux_sdata.manual_cid_name(:);
        return
    else
        % preallocate some columns for later data
        sdata = addToTable(sdata,2,'auto_cid_info',NaN(1,5),'auto_cid_type',NaN(1,1),'auto_cid_name',cell(1,1),'manual_cid',NaN(1,1),'manual_cid_name',cell(1,1)); % [table in, method 1=add 2=preallocate, variable name, variable value]    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% Otherwise we need to present a figure of each cluster to the user
%% So they can decide on the cell type
% prepare figure for manual refining
    fig_overall = figure('visible','on','Position',[50,150,1800,800]);
    fig_hor = 7; % how many plots wide should it be
    fig_ver = nsess; % how many plots tall should it be
    fspac = 0.01; % the spacing around the plots, on all sides
    fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
    fmarg = 0.03; % the margins around the plots, at the edge of the figure
    fsiz = 8; % the fontsize for different texts
    flw = 1; % the line width for different plots  
    global ctype
    ann = annotation('textbox',[0, 1, 1, 0],'string','Initializing...','FontSize',16,'LineStyle','none','interpreter','none');  

%% ##################################### Run through every cluster
    for kk = 1:length(cnames)
        % clear figure for next cluster
        clf(fig_overall);     
        delete(ann);
        uci = cnames{kk};

        % get the tetrode and cluster of the cell, skip if noise
        idx = strcmp(sdata.uci,uci);
        tet = sdata.tet(idx,1);
        clu = sdata.clu(idx,1);
        if ~clu
            continue
        end        

%% #################### automated cell typing
        big_frate = 0;
        for ss = 1:nsess
            part_now = snames{ss};
            pidx = strcmp(sdata.uci,uci) & strcmp(sdata.part,part_now);

            % get the cell's firing rate in this part
            pfrate = sdata.frate(pidx,1); % firing rate in part
            if ~pfrate || pfrate<big_frate
                continue
            end
            big_frate = pfrate;

            % get some other cell characteristics
            skaggs = sdata.spatial_info_bsec(pidx,1);
            gscore = sdata.grid_score(pidx,1);
            hscore = sdata.hd_rayleigh(pidx,1);
            wav_maxs = sdata.waveform_max(pidx,:);
            wav_wide = sdata.waveform_width(pidx,:);
            
            % find the wave width associated with highest mean amplitude
            [~,mindx] = nanmax(wav_maxs);
            max_width = wav_wide(mindx);

            cid_info = [max_width pfrate skaggs gscore hscore];
        end
        sdata.auto_cid_info(idx,:) = repmat(cid_info,sum(idx),1);
        sdata.auto_cid_type(idx,:) = 6; % preallocate with noise ID
        sdata.auto_cid_name(idx,:) = {'Noise'}; % preallocate with noise ID        
        
        cols = {'w','w','w','w','w','w','w','w'};
        if cid_info(1) > 250
            sdata.auto_cid_type(idx,:) = 2;
            sdata.auto_cid_name(idx,:) = {'Pyramidal'};
            cols{2} = 'g';            
            if cid_info(2) >0.1 && cid_info(2) <5 && cid_info(3)>0.5
                sdata.auto_cid_type(idx,:) = 1;
                sdata.auto_cid_name(idx,:) = {'Place cell'};
                cols{1} = 'g';
                cols{2} = 'w';
            end
        else
            sdata.auto_cid_type(idx,:) = 5;  
            sdata.auto_cid_name(idx,:) = {'Interneuron'}; 
            cols{5} = 'g';
        end   
    
%% #################### manual refining
        % pushbuttons for cell typing
        bxpos = 1600;
        ctype = 0;
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{1},'String','Place cell','Position',[bxpos 700 150 50],'Callback',{@change_ctype,1});  
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{2},'String','Pyramidal cell','Position',[bxpos 620 150 50],'Callback',{@change_ctype,2}); 
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{3},'String','Grid cell','Position',[bxpos 540 150 50],'Callback',{@change_ctype,3}); 
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{4},'String','Agree','Position',[bxpos 380 150 50],'Callback',{@change_ctype,-1});     
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{5},'String','Interneuron','Position',[bxpos 300 150 50],'Callback',{@change_ctype,5});    
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{6},'String','Noise','Position',[bxpos 140 150 50],'Callback',{@change_ctype,6});              
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{7},'String','Boundary cell','Position',[bxpos 220 150 50],'Callback',{@change_ctype,7});  
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{8},'String','Periodic cell','Position',[bxpos 460 150 50],'Callback',{@change_ctype,8});  

        for ss = 1:nsess
            part_now = snames{ss};
            pidx = strcmp(sdata.uci,uci) & strcmp(sdata.part,part_now);
            
            part_duration = sdata.duration(pidx,1);
            pfrate = sdata.frate(pidx,1); % firing rate in part
            ps = (fig_hor * ss) - fig_hor + 1;

            %% spike/position plot
            pox = double(pdata.(part_now).pox);
            poy = double(pdata.(part_now).poy);
            pot = double(pdata.(part_now).pot);            
            spx = pox(sdata.part_spike_index{pidx,1});
            spy = poy(sdata.part_spike_index{pidx,1});
            if numel(spx)<2 % if there are no spikes  
                continue
            end      
            
%% Spike and position plot        
            subaxis(fig_ver,fig_hor,ps,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                % plot position data, excluding pieces not included in this part
                % by inserting NaNs between intervals we can plot this as one line, which saves on memory
                dindax = abs([0; diff(pot)])>0.1;
                pos_plot = [pox poy];
                pos_plot(dindax,:) = NaN;
                plot(pos_plot(:,1),pos_plot(:,2),'Color',[.5 .5 .5]); hold on;

                % plot spikes after position data so they are all on top
                plot(spx,spy,'ro','MarkerFaceColor','r','MarkerSize',3)    

                % additional settings
                daspect([1 1 1])
                axis xy off square tight
                title(sprintf('Spk:%d t:%ds Fr:%.2f',numel(spx),round(part_duration),pfrate),'FontSize',6,'FontWeight','normal');   
                ax = gca;
                ax.XLim = [min(pox),max(pox)];
                ax.YLim = [min(poy),max(poy)];
                axis xy off  

%% Firing rate map
            subaxis(fig_ver,fig_hor,ps+1,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
                ratemap = sdata.ratemap{pidx,1};
                skaggs = sdata.spatial_info_bsec(pidx,1);
                spars = sdata.sparsity(pidx,1);
                cohe = sdata.spatial_coherence(pidx,1);
                
                im = imagesc(ratemap);
                set(im,'alphadata',~isnan(ratemap));
                daspect([1 1 1])
                axis xy off
                title(sprintf('SI: %.2f, Sp: %.2f, Co: %.2f',skaggs,spars*100,cohe),'FontSize',6,'FontWeight','normal');

%% HD plot
            subaxis(fig_ver,fig_hor,ps+2,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);  
                ai = linspace(0,2*pi,pdata.config.hd_bins)'; % angles for binning   
                hd_ratemap = sdata.hd_ratemap{pidx};
                hd_dwellmap = pdata.(part_now).hd_dwellmap;            
            
                hd_dwellmap = hd_dwellmap ./ max(hd_dwellmap) .* max(hd_ratemap);       
                mmp = mmpolar(ai,hd_dwellmap,'k',ai,hd_ratemap,'b','FontSize',fsiz,'Grid','on','RGridVisible','off','RTickVisible','off','TTickDelta',20,'RTickLabelVisible','on','TTickLabelVisible','on');
                set(mmp(1),'LineWidth',0.5);
                set(mmp(2),'LineWidth',0.5);
                p1 = patch(get(mmp(1),'XData'),get(mmp(1),'YData'),'k','FaceAlpha',0.2);
                p2 = patch(get(mmp(2),'XData'),get(mmp(2),'YData'),'b','FaceAlpha',0.5);

                ax1 = get(gca,'Position');
                ll = legend([p1 p2],{'Session','Cell'},'Units','normalized','FontSize',6,'NumColumns',2,'Location','northoutside'); legend boxoff;
                set(gca,'Position',ax1);
                legend boxoff
                axis square                
            
            
            
                
%% theta autocorrelogram
            subaxis(fig_ver,fig_hor,ps+3,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);     
                autoc = sdata.t500_spike_autocorr{pidx};
                tlag = -495:10:495;
                bar(tlag,autoc,1,'k'); hold on;
                ax = gca;
                ax.XLim = [-500 500];
                axis square
                xlabel('Time lag (ms)') % label x axis
                ylabel('Probability') % label y axis
                set(gca,'LineWidth',1,'layer','top','FontSize',fsiz); 
                
%% refractory period
            subaxis(fig_ver,fig_hor,ps+4,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);     
                t25 = histcents(-25:0.5:25);
                c25 = sdata.t25_spike_autocorr{pidx};
                bar(t25,c25,0.9,'k');
                v1 = axis;
                axis([t25(1) t25(end) v1(3) v1(4)]);
                xlabel('Time lag (ms)') % label x axis
                set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);                 
                %title(sprintf('RPV: %d',RPV),'FontSize',6);
                hold on
                v = axis;
                plot([-6; -6],[0 v(4)],'r:','LineWidth',0.5)
                plot([6; 6],[0 v(4)],'r:','LineWidth',0.5)      
                axis square
  
%% waveform
            subaxis(fig_ver,fig_hor,ps+5,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);        
                wavtime = -200:20:780;
                maxwavs = sdata.waveform_max(pidx,:);
                [maxval,maxidx] = nanmax(maxwavs);
                mwaves = reshape(cell2mat(sdata.waveform_mean(pidx,:)),1,50,4);
                swaves = reshape(cell2mat(sdata.waveform_stdv(pidx,:)),1,50,4);
                ax_lim = [floor(nanmin(mwaves(:)-swaves(:))./10).*10 ceil(nanmax(mwaves(:)+swaves(:))./10).*10];

                % plot the mean and SD of the waveform with the highest amplitude
                [hl,hp] = boundedline(wavtime,mwaves(1,:,maxidx),swaves(1,:,maxidx),'-k');
                set(hl,'Color','r','LineStyle','-','LineWidth',1); set(hp,'FaceColor','b','FaceAlpha',.5);     
                ax = gca;
                ax.XLim = [-200 780];
                ax.YLim = ax_lim;
                box on
                grid on
                axis square
                xlabel('Time (ms)')
                ylabel(sprintf('Amplitude (%cV)',char(956)))
                text(.5,1.05,sprintf('Peak: %.1f%cV, Width: %.fms, SNR: %.1f',maxval,char(956),sdata.waveform_width(pidx,maxidx),sdata.waveform_snr(pidx,maxidx)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')           
        end

%% #################### get users response and add it to sdata
        ann_str = sprintf('Automated decision was: %s - Frate: %.2f, WoW: %d, SI: %.3f',sdata.auto_cid_name{find(idx,1,'first')},sdata.auto_cid_type(find(idx,1,'first')),sdata.auto_cid_info(find(idx,1,'first'),1),sdata.auto_cid_info(find(idx,1,'first'),3));
        ann = annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',16,'LineStyle','none','interpreter','none');  
        uiwait(fig_overall); % at this point start waiting for a cell type button to be pressed

        if ctype == -1
            ctype = sdata.auto_cid_type(find(idx,1,'first'));
        end

        sdata.manual_cid_name(idx) = {'Undefined'};
        if ctype == 1
            sdata.manual_cid_name(idx) = {'Place cell'};
        elseif ctype == 2
            sdata.manual_cid_name(idx) = {'Pyramidal cell'};      
        elseif ctype == 3
            sdata.manual_cid_name(idx) = {'Grid cell'};         
        elseif ctype == 5
            sdata.manual_cid_name(idx) = {'Interneuron'};
        elseif ctype == 6
            sdata.manual_cid_name(idx) = {'Noise'};  
        elseif ctype == 7
            sdata.manual_cid_name(idx) = {'Boundary cell'};              
        elseif ctype == 8
            sdata.manual_cid_name(idx) = {'Periodic cell'};              
        end
        sdata.manual_cid(idx) = ctype;     
        ctype = 0;  

    end
    close(fig_overall)
    
    col_idx = ismember(sdata.Properties.VariableNames,{'auto_cid_name','auto_cid_type','auto_cid_info','manual_cid_name','manual_cid'});
    faux_sdata = sdata(:,col_idx);
    save(ctname,'faux_sdata');
    save(sdname,'sdata');
    save(pdname,'pdata');    
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dumb pushbutton control function
function change_ctype(~,~,value) 
    global ctype; 
    ctype = value; 
    uiresume(gcbf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A stupid subfunction for preallocating a table
function tin = addToTable(tin,meth,varargin)
    switch meth
        case {1} % we want to add info to the table
            for i=1:2:length(varargin)
                tin.(varargin{i}) = varargin{i+1}; % table in (variable name) = variable value to assign
            end
            
        case {2} % we want to preallocate column(s)
            rnum = size(tin,1);
                for i=1:2:length(varargin)
                    tin.(varargin{i}) = repmat(varargin{i+1},rnum,1); % table in (variable name) = variable value to assign
                end            
    end
end




























































