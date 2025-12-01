function [start_t,end_t] = figKEYS(mtint,pdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an mtint file and outputs key press information, it can also plot the trials if needed
%   This is a very rough prrof of concept to be used for development. Input 'kvals' is a structure to allow 
%   more flexible/varied keypress inputs for detection
%   [stimes,etimes,vals,tims] = figKEYS(mtint,disp_plot)
%
%%%%%%%% Inputs
%   mtint = an mtint structure
%   kvals = structure specifying keypress information
%       kvals.start = key value used for trial starts
%       kvals.ends = key value used for trial ends
%   skip_plot = 0 to see plots of trials, 1 to skip this
%
%%%%%%%% Outputs
%   stimes = start times of trials
%   etimes = end times of trials
%   vals = key values
%   tims = key timestamps
%
%%%%%%%% Comments
%   10/05/17 created for Eleonore's Hexamaze data
%   04/10/17 fixed error where keypresses of a different case were being missed
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
if pdata.config.trial_override
    skip_plot = 0;
end
if pdata.config.fig_key_off
    skip_plot = 1;
end

% deal with input variables
inps = {'mtint','pdata','skip_plot'};
vals = {'0','0','1'};
reqd = [1 1 0];
for ff = 1:length(inps)
    if ~exist(inps{ff},'var')
        if reqd(ff)
            error('ERROR: vargin %s missing... exiting',inps{ff});            
        end        
        eval([inps{ff} '=' vals{ff} ';']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate template
fig_trials = figure('visible','on','Position',[100 100 1400 900]);
clear wrun; global wrun; wrun = 1;

while wrun
    keyname = [pwd '\' pdata.combined_name '.txt'];
    if exist(keyname,'file') & ~pdata.config.trial_override
        xout = readtable(keyname,'FileType','text','Delimiter','tab');
    else
        % read .txt files and accumulate their data
        disp(sprintf('\t...retrieving data from files'));   
        xout = table;
        ttime = 0;
        snames = pdata.session_names;
        for kk = 1:numel(snames)      
            [count,timestamps,type,value,text_value,intervals] = saveKEY([snames{kk} '.inp'],2); % load it    
            xnow = table(count+size(xout,1),ones(size(timestamps)).*kk,repmat(snames(kk),size(timestamps)),timestamps+ttime,type,value,text_value,intervals,'VariableNames',{'count','rec','session','tstamps','type','value','text_value','intervals'});
            xout = [xout;xnow];
            ttime = max(cumsum(pdata.recording_times(1:kk)));
        end 

        % write combined text file
        xout2 = xout;
        xout2.tstamps = num2str(xout.tstamps,'%.6f');
        writetable(xout2,keyname,'FileType','text','Delimiter','tab');
    end
    pdata.config.trial_override = 0;
    
    pos = mtint.pos.xy_pixels;
    pox = pos(:,1);
    poy = -pos(:,2);
    pot = mtint.pos.ts;

    % read contents of .key file
    timestamps = xout.tstamps(:);
    type = xout.type(:);
    value = xout.value(:);

    % extract 1s and 2s
    idx = type==75; % uint32('K') = 75 which are the manual keypresses
    trial_times = timestamps(idx);
    line_index = xout.count(idx);
    key_values = value(idx); % 49s = 1s, 50s = 2s, 's' = 115, 'e' = 101

    % find start and end times
    ikeys = pdata.interval_keys;
    for kk = 1:length(ikeys)
        ikeys{kk} = uint32(num2str(ikeys{kk}));
    end

    sindx = key_values==ikeys{1};
    eindx = key_values==ikeys{2};
    disp(sprintf('\t...found %d starts and %d ends',sum(sindx),sum(eindx)));
    if sum(sindx) ~= sum(eindx)
        disp(sprintf('\t...WARNING: fix irregular keypresses'));
        skip_plot = 0;
    end

    % run through and sort out keys, also plot data if requested
    stimes = trial_times(sindx);
    start_index = line_index(sindx);
    etimes = trial_times(eindx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assessment figure
    clf(fig_trials);
    uicontrol(fig_trials,'Style','pushbutton','String','Reload','Position',[100 10 150 25],'Callback',{@change_wrun,1});  
    uicontrol(fig_trials,'Style','pushbutton','String','Quit','Position',[300 10 150 25],'Callback',{@change_wrun,0});

    if length(stimes)<6
        cnum = length(stimes);
    else
        cnum = 6;
    end
    rnum = ceil(sum(sindx) / cnum);
    start_t = NaN(size(stimes));
    end_t = NaN(size(stimes));
    for ii = 1:length(stimes) % for every start time
        snow = stimes(ii);
        enow = etimes(find(etimes > snow,1,'first'));
        if isempty(enow) % if there is no end keypress after this start, use the end of data
            enow = max(pot);
        end
        start_t(ii) = snow;
        end_t(ii) = enow;

        pindx = pot > snow & pot < enow;

        %subplot(rnum,cnum,ii)
        subaxis(rnum,cnum,ii,'Spacing',.005,'Padding',.005,'Margin',.03,'Holdaxis');
        plot(pox(pindx),poy(pindx),'k'); hold on;
        plot(pox(find(pindx,1,'first')),poy(find(pindx,1,'first')),'g.','MarkerSize',20);
        plot(pox(find(pindx,1,'last')),poy(find(pindx,1,'last')),'r.','MarkerSize',20);
      
        daspect([1 1 1]);
        axis([min(pox) max(pox) min(poy) max(poy)]);
        axis on xy
        ax = gca;
        %ax.XTick = [];
        %ax.YTick = [];

        fdnames = pdata.part_names;
        pnow = 'none';                 
        for ff = 1:length(fdnames)
            if strcmp(fdnames{ff},'FDATA') || strcmp(fdnames{ff},'a00') || strcmp(fdnames{ff},'all')
                continue
            elseif ismember(ii,pdata.part_config.Intervals{ff})
                pnow = fdnames{ff};
            end
        end

        set(gca,'LineWidth',1,'layer','top','FontSize',6); 
        title(sprintf('trial %d (%d,%s)',ii,start_index(ii),pnow),'FontWeight','normal','FontSize',9);
    end

    ann_str = sprintf('\tFound %d starts and %d ends',sum(sindx),sum(eindx));
    annotation('textbox',[0, 0.1, 1, 0],'string',ann_str,'FontSize',12,'LineStyle','none','interpreter','none');  
    
    ann_str = sprintf('Check each trial, ensure it matches the part name in subplot title, if number of starts does not match number of ends, fix this in the .txt files');
    annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',12,'LineStyle','none','interpreter','none');  
    
    ann_str = sprintf('When you have edited keyfile press ''reload'' to see the updated trials, when you are satisfied press ''quit'' to move on to klustest');
    annotation('textbox',[0, 0.07, 1, 0],'string',ann_str,'FontSize',12,'LineStyle','none','interpreter','none');  
               
    open(keyname);
    
    if skip_plot
        wrun = 0;
        break
    end

    uiwait(fig_trials); % at this point start waiting for a cell type button to be pressed
end
    close(fig_trials);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dumb pushbutton control function
function change_wrun(~,~,value) 
    global wrun; 
    wrun = value; 
    uiresume(gcbf);
end






%     if ~exist(full_keyname,'file')
%         % retrieve data
%         disp(sprintf('\t...retrieving data from mtint'));
%         pos = mtint.pos.xy_pixels;
%         pox = pos(:,1);
%         poy = pos(:,2);
%         pot = mtint.pos.ts;
%         
%         %% extract the values associated with keypresses
%         header = mtint.header;
%         dur = [0 cumsum([header.duration])];
%         xout = table;
%         for ii = 1:length(flnmroot)
%             keyname = [filepath flnmroot{ii} '.txt'];
%             xoutn = readtable(keyname);
%             if ii>1
%                 xoutn.count(:) = xoutn.count(:)+max(xout.count(:,1));        
%                 xoutn.timestamps(:) = xoutn.timestamps(:)+dur(ii);
%             end
%             xout = [xout; xoutn];
%         end
%         xout(1,:) = [];
%         xout(xout.type(:)~=75,:) = []; % uint32('K') = 75 which are the manual keypresses
%         xout.text(:) = sprintf('%c',xout.value(:));
% 
%         % save this data in a combined .key file with only the keypress data
%         writetable(xout,full_keyname,'FileType','text','Delimiter','\t');     
%     else
%         disp(sprintf('\t...retrieving data from file: %s',full_keyname));   
% 
%         xout = readtable(full_keyname,'Delimiter','\t');
%     end
% 


