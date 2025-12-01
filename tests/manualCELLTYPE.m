function sdata = manualCELLTYPE(fname,oldct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function loads the sdata stored in a data folder
%   It extracts and plots the information for each cluster in order and asks the user to specify the celkl type
%   This information is then saved into the sdata
%   function manualCELLTYPE(fname)
%
%%%%%%%% Inputs
%   fname = the directory address of the folder where manualCELLTYPE should work
%
%%%%%%%% Outputs
%   sdata will be updated to include:
%         sdata.(uci).manual_cell_type = the cell type specified byu the user, 1:5, where 1 = place cell, 2 = pyramidal cell, 3 = grid cell, 4 = border cell, 5 = Interneuron, 6 = noise or undefined
%         sdata.(uci).manual_cell_name = the text version of the above
%
%%%%%%%% Comments
%   17/01/18 created 
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('fname','var') || isempty(fname) || all(isnan(fname(:)))
    error('ERROR: manualCELLTYPE cannot function without a folder address... exiting')
end

if ~exist('oldct','var') || isempty(oldct) || all(isnan(oldct(:)))
    oldct = 1;
end

% get sdata
outname = 'kwiktint';
sdname = [fname '\klustest3\' outname '\' outname '_sdata.mat'];
sdname2 = [fname '\klustest3\' outname '\' outname '_cell_types.mat'];
load(sdname,'sdata');

% if we just want to combine previous cell type data with a new sdata
if exist(sdname2,'file') && oldct
    load(sdname2,'faux_sdata');

    f = fieldnames(faux_sdata);
    for i = 1:length(f)
        uci = f{i};
        sst = strsplit(uci,'_');
        tet = str2double(sst{3}(regexp(sst{3},'\d')));
        tetstr = ['t' num2str(tet)];    
        clu = str2double(sst{4}(regexp(sst{4},'\d')));           
        
        new_uci = ['r' sdata.rat_num '_' sdata.date '_t' num2str(tet) '_c' num2str(clu)]; 
        
        if ~isfield(sdata,new_uci)
            keyboard
        end
        if isfield(sdata,uci) && isfield(sdata,new_uci)
            sdata = rmfield(sdata,uci);
        end
        
        sdata.(new_uci).manual_cell_type = faux_sdata.(f{i}).manual_cell_type;
        sdata.(new_uci).manual_cell_name = faux_sdata.(f{i}).manual_cell_name;
    end
    
    save(sdname,'sdata');
    return
end

fig_vis = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run through clusters
disp(sprintf('Analysing cluster data...'))

% get a list of all the clusters available
start_string = ['r' sdata.rat_num '_' sdata.date]; % the start of every cluster fieldname
fdnames = fieldnames(sdata);
cindx = strfind(fdnames,start_string);
cindx(cellfun('isempty',cindx)) = {NaN}; % fill empty index values with NaN
cindx = find(cell2mat(cindx)==1); % find where the non-NaNs are
cell_names = fdnames(cindx);

fdnames = fieldnames(sdata.part_config);
cindx = strfind(fdnames,'FDATA');
cindx(cellfun('isempty',cindx)) = {0}; % fill empty index values with 0
cindx = find(cell2mat(cindx)==0); % find where the non-NaNs are
sess_names = fdnames(cindx);
nsess = length(sess_names);

fig_overall = figure('visible',fig_vis,'Position',[50,60,1400,800]);
fig_hor = 6; % how many plots wide should it be
fig_ver = nsess; % how many plots tall should it be
fspac = 0.01; % the spacing around the plots, on all sides
fpadd = 0.01; % the spacing around the plots, on all sides, this takes more space than fspac though
fmarg = 0.03; % the margins around the plots, at the edge of the figure
fsiz = 8; % the fontsize for different texts
flw = 1; % the line width for different plots  

global ctype
    
%% Run through every cluster
faux_sdata = struct;
for kk = 1:length(cell_names)
    % clear figure for next cluster
    clf(fig_overall); 
    
    % pushbuttons for cell typing
    bxpos = 1200;
    ctype = 0;
    uicontrol(fig_overall,'Style','pushbutton','String','Place cell','Position',[bxpos 700 150 50],'Callback',{@change_ctype,1});  
    uicontrol(fig_overall,'Style','pushbutton','String','Pyramidal cell','Position',[bxpos 600 150 50],'Callback',{@change_ctype,2});      
    uicontrol(fig_overall,'Style','pushbutton','String','Grid cell','Position',[bxpos 500 150 50],'Callback',{@change_ctype,3});   
    uicontrol(fig_overall,'Style','pushbutton','String','Border cell','Position',[bxpos 400 150 50],'Callback',{@change_ctype,4});   
    uicontrol(fig_overall,'Style','pushbutton','String','Interneuron','Position',[bxpos 300 150 50],'Callback',{@change_ctype,5});    
    uicontrol(fig_overall,'Style','pushbutton','String','Noise','Position',[bxpos 200 150 50],'Callback',{@change_ctype,6});     
    
    uci = cell_names{kk};
    disp(sprintf('\t...cluster %s',uci));

    % get the tetrode and cluster of the cell
    sst = strsplit(uci,'_');
    tet = str2double(sst{3}(regexp(sst{3},'\d')));
    tetstr = ['t' num2str(tet)];    
    clu = str2double(sst{4}(regexp(sst{4},'\d')));   
    
    if ~clu
        continue
    end    

    for ss = 1:nsess
        part_now = sess_names{ss};
        
        part_duration = sdata.(part_now).duration;
        pfrate = sdata.(uci).(part_now).frate;
        ps = (fig_hor * ss) - fig_hor + 1;
        
        %% spike/position plot
        subaxis(fig_ver,fig_hor,ps,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);

        pox = double(sdata.(part_now).pox);
        poy = double(sdata.(part_now).poy);
        spx = double(sdata.(uci).(part_now).spx);
        spy = double(sdata.(uci).(part_now).spy);

        if ~numel(spx) % if there are no spikes  
            continue
        end        
        
        plot(pox,poy,'k')
        hold on
        plot(spx,spy,'r.','MarkerSize',20)        
        title(sprintf('Spk:%d\nt:%ds\nFr:%.2f',numel(spx),round(part_duration),pfrate),'FontSize',6);   
        daspect([1 1 1]);
        ax = gca;
        ax.XLim = [min(pox),max(pox)];
        ax.YLim = [min(poy),max(poy)];
        axis xy off
        tt = text(min(pox)-50,min(poy)+50,part_now,'FontSize',12);
        set(tt,'rotation',90);    
        
        %% ratemap
        subaxis(fig_ver,fig_hor,ps+1,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);
        
        ratemap = double(sdata.(uci).(part_now).ratemap);
        skaggs = sdata.(uci).(part_now).spatial_measures.spatial_information;
        spars = sdata.(uci).(part_now).spatial_measures.sparsity;
        MMF = sdata.(uci).(part_now).spatial_measures.mean_method_focus;
        jsd = sdata.(uci).(part_now).spatial_measures.jensen_shannon_divergence;

        im = imagesc(ratemap);
        set(im,'alphadata',~isnan(ratemap));
        daspect([1 1 1])
        axis xy off
        set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);   
        title(sprintf('SI: %.2f, Sp: %.2f, MMF: %.2f',skaggs,spars*100,MMF),'FontSize',6);
        tt = text(-1,5,sprintf('%c: %.2f (Hz), %c: %.2f (Hz), %c: %.2f (Hz)',char(708),nanmax(ratemap(:)),char(181),nanmean(ratemap(:)),char(709),nanmin(ratemap(:))),'FontSize',6);
        set(tt,'rotation',90);        
        
        %% autocorrelation
        subaxis(fig_ver,fig_hor,ps+2,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);        
        automap = double(sdata.(uci).(part_now).grid_autocorrelation);
        if isfield(sdata.(uci).(part_now),'grid_metrics')
            grid_score = sdata.(uci).(part_now).grid_score;
            grid_wavelength = sdata.(uci).(part_now).grid_metrics.wavelength;
            grid_orientation = sdata.(uci).(part_now).grid_metrics.orientation;
        else
            grid_score = sdata.(uci).(part_now).grid_score;
            grid_wavelength = sdata.(uci).(part_now).grid_spacing;
            grid_orientation = sdata.(uci).(part_now).grid_orientation;            
        end

        im = imagesc(automap);
        set(im,'alphadata',~isnan(automap));
        daspect([1 1 1])
        axis xy off
        set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz); 
        title(sprintf('G: %.2f, Wavelength: %.2f, Orientation: %.2f',grid_score,grid_wavelength,grid_orientation),'FontSize',6);

        %% refractory period
        subaxis(fig_ver,fig_hor,ps+3,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);        
        tms2 = sdata.(uci).(part_now).refractory_period(:,2);
        corrdata2 = sdata.(uci).(part_now).refractory_period(:,1);
        RPV = sdata.(uci).(part_now).rpv_total; % add data to structure
        tau_r = sdata.tau_r;

        bar(tms2,corrdata2,0.9,'k');
        v1 = axis;
        axis([tms2(1) tms2(end) v1(3) v1(4)]);
        xlabel('Time lag (ms)') % label x axis
        set(gca,'LineWidth',flw,'layer','top','FontSize',fsiz);                 
        title(sprintf('RPV: %d',RPV),'FontSize',6);
        hold on
        v = axis;
        plot([-tau_r; -tau_r],[0 v(4)],'r','LineWidth',1)
        plot([tau_r; tau_r],[0 v(4)],'r','LineWidth',1) 
        plot([-6; -6],[0 v(4)],'r:','LineWidth',0.5)
        plot([6; 6],[0 v(4)],'r:','LineWidth',0.5)      
        axis square
        
        %% waveform
        subaxis(fig_ver,fig_hor,ps+4,'Spacing',fspac,'Padding',fpadd,'Margin',fmarg);        
        wavtime = -200:20:780;
        maxwavs = sdata.(uci).(part_now).waveform_max;
        [ax_high,mval] = nanmax(maxwavs);
        ax_high = ax_high + max(max(sdata.(uci).(part_now).waveform_stdv{mval}));
        ax_low = sdata.(uci).(part_now).waveform_min(mval) - max(max(sdata.(uci).(part_now).waveform_stdv{mval}));
        [hl,hp] = boundedline(wavtime,sdata.(uci).(part_now).waveform_mean{mval},sdata.(uci).(part_now).waveform_stdv{mval},'-k');
        set(hl,'Color','r') % line color
        set(hl,'LineStyle','-') % line style
        set(hl,'LineWidth',1) % line width
        set(hp,'FaceColor','b') % color of area
        set(hp,'FaceAlpha',1) % transparency of area
        ax = gca;
        ax.XLim = [-200 780];
        ax.YLim = [ax_low ax_high];
        box on
        xlabel('Time (ms)')
        ylabel(sprintf('Amplitude (%cV)',char(956)))
        title(sprintf('Peak: %.1f%cV, Width: %.fms, SNR: %.2f',sdata.(uci).(part_now).waveform_max(mval),char(956),sdata.(uci).(part_now).waveform_width(mval),sdata.(uci).(part_now).channel_snr(mval)),'FontSize',6);
        grid on        
        axis square
        
    end

    %% get users response and add it to sdata
    uiwait(fig_overall); % at this point start waiting for a cell type button to be pressed
    sdata.(uci).manual_cell_type = ctype;
    faux_sdata.(uci).manual_cell_type = ctype;
    
    if ctype == 1
        cname = 'Place cell';
    elseif ctype == 2
        cname = 'Pyramidal cell';        
    elseif ctype == 3
        cname = 'Grid cell';
    elseif ctype == 4
        cname = 'Border cell';
    elseif ctype == 5
        cname = 'Interneuron';
    elseif ctype == 6
        cname = 'Noise';        
    end
    sdata.(uci).manual_cell_name = cname;
    faux_sdata.(uci).manual_cell_name = ctype;    
    disp(sprintf('\t\t...is a %s',cname));
    ctype = 0;

end

disp(sprintf('\t...saving new sdata'));
save(sdname,'sdata');
save(sdname2,'faux_sdata');
close(fig_overall)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dumb pushbutton control function
function change_ctype(~,~,value) 
    global ctype; 
    ctype = value; 
    uiresume(gcbf);
end




























