function [stimes,etimes] = manageKEYS(var1,meth,interval_keys,skip_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function saves or loads a .key file (a text file version of an .inp file)
%   [count,timestamps,type,value] = saveKEY(iname,meth)
%
%%%%%%%% Inputs
%   iname = the name of the .inp file (or key file, it doesn't really matter)
%   meth = (default = 1) 1 to save a .key file for a corresponding .inp file, 2 to load the contents of an existing .key file
%
%%%%%%%% Outputs
%   count = the number of inputs
%   timestamps = the time values for these
%   type = the type of each
%   value = the value of each (uint32 integer)
%   text_value = the value of each (actual string or keyboard value)
%
%%%%%%%% Comments
%   19/05/16 created 
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if isstruct(var1) % if the first input is an mtint
    mtint = var1;
    snames = mtint.flnmroot;
elseif iscell(var1) % if the first input is a cell array of session names
    snames = var1;
    mtint = struct();
end

global open_method
open_method = 2; % set to 1 to open key files in Matlab, 2 to open them in notepad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process .klustest_keys file
% for ii = 1:length(snames)
%     [~,nme,~] = fileparts(snames{ii});
%     keyname1 = [pwd '\' nme '.klustest_keys']; 
%     keyname2 = [pwd '\' nme '.inp_keys']; 
%     if ~exist(keyname1,'file') && exist(keyname2,'file')
%         temp = readtable('Y:\Spiers lab\Eleonore\FourRoomExp\r35\2020-01-03\35_4R10a.inp_keys.txt','Delimiter','\t');
%         xin = table;
%         xin.timestamps = cellstr(num2str(temp(:,1)));
%         xin.text_value = temp(:,2);
%         writetable(xin,keyname1,'FileType','text','Delimiter','tab');        
%     end       
% end

% show a message so the user know which files will be created/overwritten
rekey = [];
revec = false(length(snames));
for ii = 1:length(snames)
    % sort out file names
    [~,nme,~] = fileparts(snames{ii});
    keyname = [pwd '\' nme '.klustest_keys']; 
    if ~exist(keyname,'file') || meth==1
        if isempty(rekey)
            rekey = [nme '.klustest_keys'];            
        else
            rekey = [rekey ', ' nme '.klustest_keys'];
        end
        revec(ii) = true;
    end       
end
disp(sprintf('\t...key files that will be created: %s',rekey))
    
% run through the data and actually make the files
for ii = 1:length(snames)
    if revec(ii) % if we need or want to (re)create the key file
        % sort out file names
        [~,nme,~] = fileparts(snames{ii});
        inpname = [pwd '\' nme '.inp'];
        keyname = [pwd '\' nme '.klustest_keys']; 

        % read .inp file            
        if ~exist(inpname,'file')
            error(sprintf('ERROR: .inp file %s not found... unable to load',inpname))
        end
        [count,timestamps,type,value] = read_key(inpname);

        % cut the data to the manual keypresses only
        mandx = uint32(cell2mat(type(:)))==75;
        if ~any(mandx)
            xin = table;
            xin.timestamps = zeros(0);
            xin.text_value = zeros(0);    
        else      
            value = value(mandx);
            timestamps = timestamps(mandx);

            % process the value data to fix errors
            value = cellfun(@transpose,value,'UniformOutput',0); 
            value = cellfun(@num2str,value,'UniformOutput',0);
            value = cellfun(@strtrim,value,'UniformOutput',0);
            value = regexprep(value,{'\W',' '},'');
            value(cellfun(@isempty,value)) = {'#'}; % some keys are not recognised by uint32, replace these with NaNs, usually this happens when we hit a weird key by accident like ` or ¬
            value = value';

            xin = table;
            xin.timestamps = cellstr(num2str(timestamps'));
            xin.text_value = cell2mat(value(:));
        end
        % write key file
        writetable(xin,keyname,'FileType','text','Delimiter','tab');
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate template
fig_trials = figure('visible','on','Position',[100 100 1400 900]);
clear wrun; global wrun; wrun = 1; global ssessions;

% snamesx = ssessions;
load_pos = 1;
while wrun    
    % load the positions, but not if we already loaded them in a previous loop
    if load_pos
        total_duration = 0;
        all_pos = [];
        all_pot = [];
        if ~isempty(fieldnames(mtint)) % if we were given an mtint, use the positions there to save time
            disp(sprintf('\t...loading positions from mtint'))
            all_pos = mtint.pos.led_pos;
            all_pot = mtint.pos.ts;      

        else % if we did not get an mtint
            disp(sprintf('\t...loading positions from .pos files'))        
            loopout = looper(length(snames));
            for ii = 1:length(snames)      
                h2 = getDACQDATAHEADERS(snames(ii),'.set');
                [pos,post] = rawpos([snames{ii},'.pos'],h2.colactive_1+h2.colactive_2); % raw led positions
                all_pos = [all_pos; pos];
                all_pot = [all_pot; post(:)+total_duration];
                h1 = getDACQDATAHEADERS(snames(ii),'.pos');
                total_duration = total_duration + h1.duration;
                loopout = looper(loopout);
            end 
        end        
        load_pos = 0;
    end
    
    % load the keypresses
    disp(sprintf('\t...loading keypresses'))
    xout = table;    
    loopout = looper(length(snames));
    total_duration = 0;
    for ii = 1:length(snames)      
        [~,nme,~] = fileparts(snames{ii});
        keyname = [pwd '\' nme '.klustest_keys'];     
        xnow = readtable(keyname,'FileType','text','Delimiter','tab');
        if ~isempty(xnow)
            xnow.timestamps = xnow.timestamps+total_duration;
            xnow.timestamps2 = xnow.timestamps;        
            xnow.session(:) = snames(ii);
            xout = [xout; xnow];
        end
        h1 = getDACQDATAHEADERS(snames(ii),'.pos');
        total_duration = total_duration + h1.duration; 
        loopout = looper(loopout);
    end         

    % read contents of .key file
    trial_times = xout.timestamps;
    key_values = strtrim(xout.text_value);

    % find start and end times
    ikeys = interval_keys;

    sindx = ismember(key_values,ikeys{1});
    eindx = ismember(key_values,ikeys{2});
    disp(sprintf('\t...found %d starts and %d ends',sum(sindx),sum(eindx)));
    bcol = 'g';
    if sum(sindx) ~= sum(eindx)
        disp(sprintf('\t...WARNING: fix irregular keypresses'));
        skip_plot = 0;
        bcol = 'r';
    end

    % run through and sort out keys, also plot data if requested
    stimes = trial_times(sindx);
    stimes2 = xout.timestamps2;
    stimes2 = stimes2(sindx);    
    ssessions = {xout.session{sindx}}';
    etimes = trial_times(eindx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assessment figure
    clf(fig_trials);
    uicontrol(fig_trials,'Style','pushbutton','BackgroundColor','g','String','Reload','Position',[100 10 150 25],'Callback',{@change_wrun,1});  
    uicontrol(fig_trials,'Style','pushbutton','BackgroundColor',bcol,'String','Quit','Position',[300 10 150 25],'Callback',{@change_wrun,0});

    if length(stimes)<6
        cnum = length(stimes);
    else
        cnum = 6;
    end
    rnum = max([ceil(sum(sindx) / cnum) 2]);
    start_t = NaN(size(stimes));
    end_t = NaN(size(stimes));
    for ii = 1:length(stimes) % for every start time
        snow = stimes(ii);
        enow = etimes(find(etimes > snow,1,'first'));
        if isempty(enow) % if there is no end keypress after this start, use the end of data
            enow = max(all_pot);
        end
        start_t(ii) = snow;
        end_t(ii) = enow;

        pindx = all_pot > snow & all_pot < enow;

        subaxis(rnum,cnum,ii,'Spacing',.005,'Padding',.005,'Margin',.03,'Holdaxis');
        pox = all_pos(:,1,1);
        poy = all_pos(:,1,2);
        plot(pox(pindx),poy(pindx),'k'); hold on;
        plot(pox(find(pindx,1,'first')),poy(find(pindx,1,'first')),'g.','MarkerSize',20);
        plot(pox(find(pindx,1,'last')),poy(find(pindx,1,'last')),'r.','MarkerSize',20);
      
        daspect([1 1 1]);
        axis([min(pox) max(pox) min(poy) max(poy)]);
        axis on xy
        ax = gca;
        %ax.XTick = [];
        %ax.YTick = [];
        title(sprintf('Session: %s\nStart t: %.2fs',ssessions{ii},stimes2(ii)),'FontWeight','normal','HorizontalAlignment','left','FontSize',6);

        set(gca,'units','pixels')
        set(gcf,'units','pixels')
        ax = gca;
        uicontrol(gcf,'Style','pushbutton','String','Open','Units','pixels','Position',[ax.Position(1) ax.Position(2) 50 25],'Callback',{@open_file,ii});  
        set(gca,'LineWidth',1,'layer','top','FontSize',6); 
    end

    ann_str = sprintf('\tFound %d starts and %d ends',sum(sindx),sum(eindx));
    annotation('textbox',[0, 0.1, 1, 0],'string',ann_str,'FontSize',12,'LineStyle','none','interpreter','none');  
    
    ann_str = sprintf('Check each trial, ensure it matches the part name in subplot title, if number of starts does not match number of ends, fix this in the .txt files AND then reload');
    annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',12,'LineStyle','none','interpreter','none');  
    
    ann_str = sprintf('When you have edited keyfile press ''reload'' to see the updated trials, when you are satisfied press ''quit'' to move on to klustest');
    annotation('textbox',[0, 0.07, 1, 0],'string',ann_str,'FontSize',12,'LineStyle','none','interpreter','none');  
    
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

function open_file(~,~,value) 
    global ssessions
    global open_method
    if open_method==1 % open in Matlab
        open([pwd '\' ssessions{value} '.klustest_keys']);
    elseif open_method==2 % open in notepad
        system(['notepad', ' ', fullfile(pwd,[ssessions{value},'.klustest_keys', '&'])]);
    end
end




































    










