
fig_vis = 'off'; % whether you want figure to be visible or not




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get information about the parts and intervals
% [fname,pname] = uigetfile('*.mat','Select sdata.mat');

disp(sprintf('Loading sdata...'))
% load([pname fname]);

% part_config is saved inside sdata, we need it to get the interval times etc
disp(sprintf('Extracting part_config...'))
part_config = sdata.part_config;

% find which parts were separated by intervals (i.e. part_config method = 3)
pnames = fieldnames(part_config);
mnames = {};
mtimes = {};
for pp = 1:length(pnames) % for every part
    pname_now = pnames{pp,1};
    
    % if this is the FDATA field just skip it
    if strcmp(pname_now,'FDATA')
        continue
    end
    
    % if the part's method is 3, add it to mnames cell array and add its time values to mtimes cell array
    meth_now = part_config.(pname_now).method;
    if meth_now == 3
        mnames{pp,1} = pname_now; % part name in structure
        mtimes{pp,1} = part_config.(pname_now).times; % interval times
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now we can go and get cell data
% every cell is saved in a field which has a name starting with: r (rat num) _ date
% so we can find fields with this in their name to get all cells
start_string = ['r' sdata.rat_num '_' sdata.date];
anames = fieldnames(sdata);
cindx = strfind(anames,start_string);
cindx(cellfun('isempty',cindx)) = {0};
cindx = logical(cell2mat(cindx));
cell_names = anames(cindx); % these are all the cell names

for cc = 1:length(cell_names) % for every cluster/cell
    cnow = cell_names{cc,1}; % the current cell name
    
    % if this is a noise cluster just skip it
    if strfind(cnow,'c0')
        continue
    end
    
    for pp = 1:length(mnames) % for every part that was separated by method 3 (intervals)
        mnow = mnames{pp,1}; % the current part name
        mtime = mtimes{pp,1}; % the current part interval times
        
        % get position data for this part
        pox = sdata.(mnow).pox; % position x coordinates
        poy = sdata.(mnow).poy; % position x coordinates
        pot = sdata.(mnow).pot; % position x coordinates

        % get the spike data for this cell
        spx = sdata.(cnow).(mnow).spx; % spike x coordinates
        spy = sdata.(cnow).(mnow).spy; % spike y coordinates
        spt = sdata.(cnow).(mnow).spt; % spike times
        
        % create a figure to display data
        fig_intervals = figure('visible',fig_vis,'Position',[100, 100, 1200, 800]);
        set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
        set(gcf,'color','w'); % makes the background colour white
        colormap(jet(256)); % to make sure the colormap is not the horrible default one
        fspac = 0.001; % the spacing around the plots, on all sides
        fpadd = 0.001; % the spacing around the plots, on all sides, this takes more space than fspac though
        fmarg = 0.03; % the margins around the plots, at the edge of the figure
        fsiz = 6; % the fontsize for different texts
        flw = 1; % the line width for different plots
        fig_hor = 10; % how many plots wide should it be, change this if you like
        fig_ver = ceil(length(mtime(:,1))/fig_hor); % how many plots tall should it be
        
        % add an annotation to the figure with some important info
        ann_str = sprintf('Cell: %s, Part: %s, Analysed: %s',cnow,mnow,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',fsiz,'LineStyle','none','interpreter','none');      
        
        for ii = 1:length(mtime(:,1)) % for every interval
            % cut the spike and position data to the current interval time
            i_time = mtime(ii,:); % get the start and end time
            sindx = find(spt > i_time(1) & spt < i_time(2)); % find the spikes falling into this interval
            pindx = find(pot > i_time(1) & pot < i_time(2)); % find the position data falling into this interval

            poxi = pox(pindx);
            poyi = poy(pindx);
            poti = pot(pindx);
            
            spxi = spx(sindx);
            spyi = spy(sindx);
            spti = spt(sindx);

            subaxis(fig_ver,fig_hor,ii,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg,'Holdaxis');
            hold on
            plot(poxi,poyi,'k');
            plot(spxi,spyi,'r.','MarkerSize',15);
            daspect([1 1 1])
            axis([min(pox) max(pox) min(poy) max(poy)]);
            title(sprintf('I%d %.f-%.f',ii,i_time(1),i_time(2)),'FontSize',fsiz)
            axis on
            ax = gca;
            ax.XTickLabel = [];
            ax.YTickLabel = [];  
        end
        
        figfile = [pwd '\intervals\'];
        [~,~,~] = mkdir(figfile);
        print(fig_intervals,'-dpng','-r150',[figfile cnow '_' mnow '.png'])
        
        
        
    end
    

    
end



























