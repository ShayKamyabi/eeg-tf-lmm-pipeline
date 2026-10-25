clear;
close all;
clc;

tf_file = 'F:\my data\Figures\TF_MASTER_SIMPLE_FINAL\TF_OUTPUT_FINAL.mat';
save_dir = 'F:\my data\Figures\TF_MASTER_SIMPLE_FINAL\Publication_FINAL_v2';

if ~exist(tf_file,'file')
    error(['Canonical TF file not found:\n%s\n\n' ...
        'Edit tf_file at the top of this script.'],tf_file);
end

if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

HIGH_IDS = [ ...
     6 17 21 22 25 29 31 32 33 36 ...
    37 39 41 46 49 54 55 58 61];

LOW_IDS = [ ...
    1  2 11 12 13 15 27 28 38 40 ...
    45 47 48 50 52 53 60 66 67];

beta_range = [13 30];
analysis_window = [0 1000];
plot_window = [-200 1200];

cascade_windows = [ ...
    200 400;
    400 600;
    600 800;
    800 1000];

cascade_labels = { ...
    '200-400 ms', ...
    '400-600 ms', ...
    '600-800 ms', ...
    '800-1000 ms'};

roi_labels = {'F3','F4','Fz','FCz','P3','P4','Pz','CPz','POz'};
exclude_topo_channels = {'A1','A2'};

% Light smoothing for display only.
% Source data and statistics remain unchanged.
display_smooth_kernel = [1 2 1] ./ 4;   % very light smoothing
display_interp_factor = 8;              % denser curve for smoother display

COL_RP   = [0 114 178]/255;
COL_HA   = [213 94 0]/255;
COL_LOW  = [60 60 60]/255;
COL_HIGH = [0 114 178]/255;
COL_ZERO = [0.60 0.60 0.60];
COL_AXIS = [0.35 0.35 0.35];
COL_WINDOW = [0.95 0.95 0.95];

COMMON_CMAP = jet(256);

fprintf('Loading TF data\n');
fprintf('%s\n',tf_file);

S = load(tf_file);

requiredFields = { ...
    'grand_ersp_RP', ...
    'grand_ersp_HA', ...
    'times', ...
    'freqs', ...
    'channels', ...
    'chanlocs_selected', ...
    'subject_id'};

for k = 1:numel(requiredFields)
    if ~isfield(S,requiredFields{k})
        error('Missing required field "%s" in TF_OUTPUT_FINAL.mat.', requiredFields{k});
    end
end

RP = S.grand_ersp_RP;
HA = S.grand_ersp_HA;
times = double(S.times(:));
freqs = double(S.freqs(:));
chanNames = cellstr(string(S.channels(:)));
subject_id = double(S.subject_id(:));
chanlocs_saved = S.chanlocs_selected;

if isfield(S,'subject_names')
    subject_names = string(S.subject_names(:));
else
    subject_names = strings(numel(subject_id),1);
end

if isfield(S,'have_chanloc')
    have_chanloc_saved = logical(S.have_chanloc(:));
else
    have_chanloc_saved = true(numel(chanNames),1);
end

N = numel(subject_id);
nCh = numel(chanNames);
nF = numel(freqs);
nT = numel(times);

if ~isequal(size(RP),[nF nT nCh N])
    error('RP dimensions do not match [freq time channel subject].');
end
if ~isequal(size(HA),[nF nT nCh N])
    error('HA dimensions do not match [freq time channel subject].');
end

fprintf('TF dimensions = %d freq x %d time x %d channels x %d subjects\n', nF,nT,nCh,N);
fprintf('Time range = %.3f to %.3f ms\n',min(times),max(times));
fprintf('Frequency range = %.3f to %.3f Hz\n',min(freqs),max(freqs));

idxHigh = ismember(subject_id,HIGH_IDS);
idxLow  = ismember(subject_id,LOW_IDS);

fprintf('Group check\n');
fprintf('High CC IDs found = %d / 19\n',sum(idxHigh));
fprintf('Low CC IDs found  = %d / 19\n',sum(idxLow));

fMask = freqs >= beta_range(1) & freqs <= beta_range(2);
tAnalysis = times >= analysis_window(1) & times <= analysis_window(2);
tPlot = times >= plot_window(1) & times <= plot_window(2);

if ~any(fMask)
    error('No frequency bins found in 13-30 Hz.');
end
if ~any(tAnalysis)
    error('No time bins found in 0-1000 ms.');
end
if ~any(tPlot)
    error('No time bins found in requested plotting window.');
end

fprintf('TF bins\n');
fprintf('Beta bins = %d | actual %.3f to %.3f Hz\n', ...
    sum(fMask),freqs(find(fMask,1,'first')),freqs(find(fMask,1,'last')));
fprintf('0-1000 ms bins = %d | actual %.3f to %.3f ms\n', ...
    sum(tAnalysis),times(find(tAnalysis,1,'first')),times(find(tAnalysis,1,'last')));

roi_idx = NaN(1,numel(roi_labels));
for r = 1:numel(roi_labels)
    hit = find(strcmpi(chanNames,roi_labels{r}));
    if numel(hit) ~= 1
        error('ROI9 channel "%s" matched %d channels; expected exactly 1.', roi_labels{r},numel(hit));
    end
    roi_idx(r) = hit;
end

if numel(chanlocs_saved) ~= nCh
    error('chanlocs_selected length does not match channel dimension.');
end

savedLocLabels = string({chanlocs_saved.labels});
chanlocs_ordered = chanlocs_saved([]);
haveLoc = false(nCh,1);

for c = 1:nCh
    hit = find(strcmpi(savedLocLabels,chanNames{c}));
    if numel(hit) ~= 1
        error('Could not uniquely match chanloc for channel %s.',chanNames{c});
    end
    chanlocs_ordered(c) = chanlocs_saved(hit);

    if numel(have_chanloc_saved) == nCh
        haveLoc(c) = have_chanloc_saved(c);
    else
        haveLoc(c) = true;
    end
end

topoBaseKeep = haveLoc;
for k = 1:numel(exclude_topo_channels)
    topoBaseKeep = topoBaseKeep & ~strcmpi(chanNames(:),exclude_topo_channels{k});
end

completeRP9 = false(N,1);
completeHA9 = false(N,1);
for s = 1:N
    okRP = true;
    okHA = true;
    for r = 1:numel(roi_idx)
        xRP = RP(fMask,tAnalysis,roi_idx(r),s);
        xHA = HA(fMask,tAnalysis,roi_idx(r),s);
        if ~any(isfinite(xRP(:)))
            okRP = false;
        end
        if ~any(isfinite(xHA(:)))
            okHA = false;
        end
    end
    completeRP9(s) = okRP;
    completeHA9(s) = okHA;
end
completePAIR9 = completeRP9 & completeHA9;

fprintf('ROI9 completeness\n');
fprintf('RP complete ROI9      = %d/%d\n',sum(completeRP9),N);
fprintf('HA complete ROI9      = %d/%d\n',sum(completeHA9),N);
fprintf('RP+HA paired complete = %d/%d\n',sum(completePAIR9),N);

RP_roi = NaN(nT,N);
HA_roi = NaN(nT,N);

for s = 1:N
    if completeRP9(s)
        for t = 1:nT
            x = RP(fMask,t,roi_idx,s);
            xf = x(isfinite(x));
            if ~isempty(xf)
                RP_roi(t,s) = mean(xf);
            end
        end
    end
    if completeHA9(s)
        for t = 1:nT
            x = HA(fMask,t,roi_idx,s);
            xf = x(isfinite(x));
            if ~isempty(xf)
                HA_roi(t,s) = mean(xf);
            end
        end
    end
end

DiffROI = HA_roi - RP_roi;

[meanRP,semRP,nRP_t] = local_mean_sem(RP_roi(:,completeRP9));
[meanHA,semHA,nHA_t] = local_mean_sem(HA_roi(:,completeHA9));

lowPairMask = idxLow & completePAIR9;
highPairMask = idxHigh & completePAIR9;

[meanLowDiff,semLowDiff,nLow_t] = local_mean_sem(DiffROI(:,lowPairMask));
[meanHighDiff,semHighDiff,nHigh_t] = local_mean_sem(DiffROI(:,highPairMask));

RP_scalp = NaN(nCh,N);
HA_scalp = NaN(nCh,N);

for s = 1:N
    for c = 1:nCh
        x = RP(fMask,tAnalysis,c,s);
        xf = x(isfinite(x));
        if ~isempty(xf)
            RP_scalp(c,s) = mean(xf);
        end

        x = HA(fMask,tAnalysis,c,s);
        xf = x(isfinite(x));
        if ~isempty(xf)
            HA_scalp(c,s) = mean(xf);
        end
    end
end

Diff_scalp = HA_scalp - RP_scalp;

[mapLowRP,nLowRP_ch] = local_channel_mean(RP_scalp(:,idxLow));
[mapLowHA,nLowHA_ch] = local_channel_mean(HA_scalp(:,idxLow));
[mapHighRP,nHighRP_ch] = local_channel_mean(RP_scalp(:,idxHigh));
[mapHighHA,nHighHA_ch] = local_channel_mean(HA_scalp(:,idxHigh));
[mapGrandDiff,nGrandDiff_ch] = local_channel_mean(Diff_scalp);

nWindows = size(cascade_windows,1);
CascadeLow = NaN(nCh,nWindows);
CascadeHigh = NaN(nCh,nWindows);
CascadeGrand = NaN(nCh,nWindows);
CascadeLowN = zeros(nCh,nWindows);
CascadeHighN = zeros(nCh,nWindows);
CascadeGrandN = zeros(nCh,nWindows);
windowCounts = zeros(nWindows,1);
windowActual = NaN(nWindows,2);

for w = 1:nWindows
    tStart = cascade_windows(w,1);
    tEnd = cascade_windows(w,2);

    if w < nWindows
        tMaskW = times >= tStart & times < tEnd;
    else
        tMaskW = times >= tStart & times <= tEnd;
    end

    windowCounts(w) = sum(tMaskW);
    windowActual(w,:) = [times(find(tMaskW,1,'first')) times(find(tMaskW,1,'last'))];

    windowDiff = NaN(nCh,N);
    for s = 1:N
        for c = 1:nCh
            xRP = RP(fMask,tMaskW,c,s);
            xHA = HA(fMask,tMaskW,c,s);
            rpF = xRP(isfinite(xRP));
            haF = xHA(isfinite(xHA));
            if ~isempty(rpF) && ~isempty(haF)
                windowDiff(c,s) = mean(haF) - mean(rpF);
            end
        end
    end

    [CascadeLow(:,w),CascadeLowN(:,w)] = local_channel_mean(windowDiff(:,idxLow));
    [CascadeHigh(:,w),CascadeHighN(:,w)] = local_channel_mean(windowDiff(:,idxHigh));
    [CascadeGrand(:,w),CascadeGrandN(:,w)] = local_channel_mean(windowDiff);
end

allMapMatrix = [ ...
    mapLowRP, mapLowHA, mapHighRP, mapHighHA, ...
    mapGrandDiff, CascadeGrand, CascadeLow, CascadeHigh];

finiteEveryMap = all(isfinite(allMapMatrix),2);
topoKeep = topoBaseKeep & finiteEveryMap;

if sum(topoKeep) < 20
    error('Too few common finite channels remain for reliable topoplot display.');
end

topoNames = chanNames(topoKeep);
topoLocs = chanlocs_ordered(topoKeep);

mapLowRP_topo = mapLowRP(topoKeep);
mapLowHA_topo = mapLowHA(topoKeep);
mapHighRP_topo = mapHighRP(topoKeep);
mapHighHA_topo = mapHighHA(topoKeep);
mapGrandDiff_topo = mapGrandDiff(topoKeep);
CascadeLow_topo = CascadeLow(topoKeep,:);
CascadeHigh_topo = CascadeHigh(topoKeep,:);
CascadeGrand_topo = CascadeGrand(topoKeep,:);

roiTopoIdx = NaN(1,numel(roi_labels));
for r = 1:numel(roi_labels)
    hit = find(strcmpi(topoNames,roi_labels{r}));
    if numel(hit) ~= 1
        error('ROI channel %s is unavailable in final common topo channel set.',roi_labels{r});
    end
    roiTopoIdx(r) = hit;
end

ALL_TOPO_VALUES = [ ...
    mapLowRP_topo(:); mapLowHA_topo(:); mapHighRP_topo(:); mapHighHA_topo(:); ...
    mapGrandDiff_topo(:); CascadeGrand_topo(:); CascadeLow_topo(:); CascadeHigh_topo(:)];
ALL_TOPO_VALUES = ALL_TOPO_VALUES(isfinite(ALL_TOPO_VALUES));
COMMON_LIMIT = max(abs(ALL_TOPO_VALUES));
COMMON_LIMIT = ceil(COMMON_LIMIT*10)/10;
COMMON_CLIM = [-COMMON_LIMIT COMMON_LIMIT];

fprintf('Topography color limits\n');
fprintf('Colormap = JET(256)\n');
fprintf('COMMON_CLIM = [%.2f %.2f] dB\n',COMMON_CLIM(1),COMMON_CLIM(2));

x_raw = times(tPlot);

rp_m_raw = meanRP(tPlot);   rp_s_raw = semRP(tPlot);
ha_m_raw = meanHA(tPlot);   ha_s_raw = semHA(tPlot);
ld_m_raw = meanLowDiff(tPlot);  ld_s_raw = semLowDiff(tPlot);
hd_m_raw = meanHighDiff(tPlot); hd_s_raw = semHighDiff(tPlot);

[x_plot, rp_m_plot, rp_s_plot] = local_soft_curve(x_raw,rp_m_raw,rp_s_raw,display_smooth_kernel,display_interp_factor);
[~,      ha_m_plot, ha_s_plot] = local_soft_curve(x_raw,ha_m_raw,ha_s_raw,display_smooth_kernel,display_interp_factor);
[~,      ld_m_plot, ld_s_plot] = local_soft_curve(x_raw,ld_m_raw,ld_s_raw,display_smooth_kernel,display_interp_factor);
[~,      hd_m_plot, hd_s_plot] = local_soft_curve(x_raw,hd_m_raw,hd_s_raw,display_smooth_kernel,display_interp_factor);

fig1 = figure('Color','w','Units','centimeters','Position',[0.5 0.5 31 17],'Renderer','painters');
tl1 = tiledlayout(fig1,2,5,'TileSpacing','compact','Padding','compact');
% Leave room for the colorbar.
tl1.OuterPosition = [0.00 0.00 0.91 1.00];

% RP vs HA
ax1 = nexttile(tl1,[1 2]);
hold(ax1,'on');
local_fill_sem(ax1,x_plot,rp_m_plot,rp_s_plot,COL_RP,0.14);
local_fill_sem(ax1,x_plot,ha_m_plot,ha_s_plot,COL_HA,0.14);
plot(ax1,x_plot,rp_m_plot,'Color',COL_RP,'LineWidth',2.4);
plot(ax1,x_plot,ha_m_plot,'Color',COL_HA,'LineWidth',2.4);
yl = ylim(ax1);
hwin1 = patch(ax1,[analysis_window(1) analysis_window(2) analysis_window(2) analysis_window(1)], ...
    [yl(1) yl(1) yl(2) yl(2)],COL_WINDOW, ...
    'EdgeColor','none','FaceAlpha',0.10,'HandleVisibility','off');
uistack(hwin1,'bottom');
plot(ax1,x_plot,rp_m_plot,'Color',COL_RP,'LineWidth',2.4);
plot(ax1,x_plot,ha_m_plot,'Color',COL_HA,'LineWidth',2.4);
xline(ax1,0,'-','Color',COL_AXIS,'LineWidth',0.9);
yline(ax1,0,'--','Color',COL_ZERO,'LineWidth',0.9);
xlim(ax1,[max(plot_window(1),min(times)) min(plot_window(2),max(times))]);
xlabel(ax1,'Time from stimulus onset (ms)');
ylabel(ax1,'Beta power (dB)');
title(ax1,'Frontoparietal beta dynamics','FontWeight','bold');
legend(ax1,{'Risky Positive','High Ambiguity'},'Location','best','Box','off','FontSize',8);
set(ax1,'FontName','Arial','FontSize',9,'LineWidth',0.9,'TickDir','out','Box','off','Layer','top');
text(ax1,-0.10,1.05,'A','Units','normalized','FontName','Arial','FontSize',12,'FontWeight','bold');

% Low vs High HA-RP
ax2 = nexttile(tl1,[1 2]);
hold(ax2,'on');
local_fill_sem(ax2,x_plot,ld_m_plot,ld_s_plot,COL_LOW,0.14);
local_fill_sem(ax2,x_plot,hd_m_plot,hd_s_plot,COL_HIGH,0.14);
plot(ax2,x_plot,ld_m_plot,'Color',COL_LOW,'LineWidth',2.4);
plot(ax2,x_plot,hd_m_plot,'Color',COL_HIGH,'LineWidth',2.4);
yl = ylim(ax2);
hwin2 = patch(ax2,[analysis_window(1) analysis_window(2) analysis_window(2) analysis_window(1)], ...
    [yl(1) yl(1) yl(2) yl(2)],COL_WINDOW, ...
    'EdgeColor','none','FaceAlpha',0.10,'HandleVisibility','off');
uistack(hwin2,'bottom');
plot(ax2,x_plot,ld_m_plot,'Color',COL_LOW,'LineWidth',2.4);
plot(ax2,x_plot,hd_m_plot,'Color',COL_HIGH,'LineWidth',2.4);
xline(ax2,0,'-','Color',COL_AXIS,'LineWidth',0.9);
yline(ax2,0,'--','Color',COL_ZERO,'LineWidth',0.9);
xlim(ax2,[max(plot_window(1),min(times)) min(plot_window(2),max(times))]);
xlabel(ax2,'Time from stimulus onset (ms)');
ylabel(ax2,'HA - RP beta power (dB)');
title(ax2,'Ambiguity-related beta modulation','FontWeight','bold');
legend(ax2,{'Low CC','High CC'},'Location','best','Box','off','FontSize',8);
set(ax2,'FontName','Arial','FontSize',9,'LineWidth',0.9,'TickDir','out','Box','off','Layer','top');
text(ax2,-0.10,1.05,'B','Units','normalized','FontName','Arial','FontSize',12,'FontWeight','bold');

% Unused tile
axBlank = nexttile(tl1,5);
axis(axBlank,'off');

% Grand HA-RP
axC = nexttile(tl1,6);
local_draw_topo(axC,mapGrandDiff_topo,topoLocs,roiTopoIdx,COMMON_CLIM,COMMON_CMAP,true);
title(axC,{'0-1000 ms','Grand HA - RP'},'FontSize',9.5,'FontWeight','bold');
text(axC,-0.12,1.08,'C','Units','normalized','FontName','Arial','FontSize',12,'FontWeight','bold');

% Grand cascade
for w = 1:nWindows
    ax = nexttile(tl1,6+w);
    local_draw_topo(ax,CascadeGrand_topo(:,w),topoLocs,[],COMMON_CLIM,COMMON_CMAP,false);
    title(ax,cascade_labels{w},'FontSize',9.5,'FontWeight','bold');
    text(ax,-0.12,1.08,char('C'+w),'Units','normalized','FontName','Arial','FontSize',12,'FontWeight','bold');
end

local_fixed_colorbar( ...
    fig1,COMMON_CMAP,COMMON_CLIM, ...
    [0.947 0.135 0.0085 0.255], ...
    'Beta power (dB)',8.5,9.5);
colormap(fig1,COMMON_CMAP);

try
    sgtitle(fig1,'Beta dynamics and scalp distribution of ambiguity-related modulation', ...
        'FontName','Arial','FontSize',12,'FontWeight','normal');
catch
end

local_export_figure(fig1,save_dir,'Beta_FINAL_01_Waveform_Topography_Summary_v2');

fig2 = figure('Color','w','Units','centimeters','Position',[1 1 22 15],'Renderer','painters');
tl2 = tiledlayout(fig2,2,2,'TileSpacing','compact','Padding','compact');
% Leave room for the colorbar.
tl2.OuterPosition = [0.00 0.00 0.90 1.00];
title(tl2,{'Absolute scalp beta power','13-30 Hz, 0-1000 ms'},'FontName','Arial','FontSize',12,'FontWeight','normal');

ax21 = nexttile(tl2,1);
local_draw_topo(ax21,mapLowRP_topo,topoLocs,roiTopoIdx,COMMON_CLIM,COMMON_CMAP,true);
title(ax21,{'Low CC','Risky Positive'},'FontSize',10.5,'FontWeight','bold');
text(ax21,-0.10,1.08,'A','Units','normalized','FontName','Arial','FontSize',12,'FontWeight','bold');

ax22 = nexttile(tl2,2);
local_draw_topo(ax22,mapLowHA_topo,topoLocs,roiTopoIdx,COMMON_CLIM,COMMON_CMAP,true);
title(ax22,{'Low CC','High Ambiguity'},'FontSize',10.5,'FontWeight','bold');
text(ax22,-0.10,1.08,'B','Units','normalized','FontName','Arial','FontSize',12,'FontWeight','bold');

ax23 = nexttile(tl2,3);
local_draw_topo(ax23,mapHighRP_topo,topoLocs,roiTopoIdx,COMMON_CLIM,COMMON_CMAP,true);
title(ax23,{'High CC','Risky Positive'},'FontSize',10.5,'FontWeight','bold');
text(ax23,-0.10,1.08,'C','Units','normalized','FontName','Arial','FontSize',12,'FontWeight','bold');

ax24 = nexttile(tl2,4);
local_draw_topo(ax24,mapHighHA_topo,topoLocs,roiTopoIdx,COMMON_CLIM,COMMON_CMAP,true);
title(ax24,{'High CC','High Ambiguity'},'FontSize',10.5,'FontWeight','bold');
text(ax24,-0.10,1.08,'D','Units','normalized','FontName','Arial','FontSize',12,'FontWeight','bold');

local_fixed_colorbar( ...
    fig2,COMMON_CMAP,COMMON_CLIM, ...
    [0.935 0.165 0.009 0.66], ...
    'Beta power (dB)',9,10);
colormap(fig2,COMMON_CMAP);

local_export_figure(fig2,save_dir,'Beta_FINAL_02_AbsoluteTopographies_v2');

fig3 = figure('Color','w','Units','centimeters','Position',[0.5 1 27 13.5],'Renderer','painters');
tl3 = tiledlayout(fig3,2,4,'TileSpacing','compact','Padding','compact');
% Leave room for the colorbar.
tl3.OuterPosition = [0.00 0.00 0.91 1.00];
title(tl3,{ ...
    'Spatiotemporal ambiguity-related beta modulation'; ...
    'Every scalp map = High Ambiguity - Risky Positive (HA - RP)'}, ...
    'FontName','Arial','FontSize',12,'FontWeight','normal');

for w = 1:nWindows
    ax = nexttile(tl3,w);
    local_draw_topo(ax,CascadeLow_topo(:,w),topoLocs,[],COMMON_CLIM,COMMON_CMAP,false);
    title(ax,cascade_labels{w},'FontSize',10.2,'FontWeight','bold');
    if w == 1
        text(ax,-0.25,0.50,{'Low CC';'HA - RP'},'Units','normalized','Rotation',90, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontName','Arial','FontSize',10,'FontWeight','bold','Clipping','off');
    end
end

for w = 1:nWindows
    ax = nexttile(tl3,4+w);
    local_draw_topo(ax,CascadeHigh_topo(:,w),topoLocs,[],COMMON_CLIM,COMMON_CMAP,false);
    if w == 1
        text(ax,-0.25,0.50,{'High CC';'HA - RP'},'Units','normalized','Rotation',90, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontName','Arial','FontSize',10,'FontWeight','bold','Clipping','off');
    end
end

local_fixed_colorbar( ...
    fig3,COMMON_CMAP,COMMON_CLIM, ...
    [0.940 0.180 0.0085 0.620], ...
    'HA - RP beta power (dB)',9,10);
colormap(fig3,COMMON_CMAP);

local_export_figure(fig3,save_dir,'Beta_FINAL_03_SpatiotemporalCascade_v2');

TimecourseTable = table( ...
    times, meanRP,semRP,nRP_t, meanHA,semHA,nHA_t, meanLowDiff,semLowDiff,nLow_t, meanHighDiff,semHighDiff,nHigh_t, ...
    'VariableNames',{'Time_ms','RP_Mean','RP_SEM','RP_N','HA_Mean','HA_SEM','HA_N', ...
    'LowCC_HAminusRP_Mean','LowCC_HAminusRP_SEM','LowCC_HAminusRP_N', ...
    'HighCC_HAminusRP_Mean','HighCC_HAminusRP_SEM','HighCC_HAminusRP_N'});
writetable(TimecourseTable,fullfile(save_dir,'Beta_FINAL_Timecourses_SourceData_v2.csv'));

AbsoluteTopoTable = table( ...
    string(chanNames(:)), mapLowRP,mapLowHA,mapHighRP,mapHighHA, nLowRP_ch,nLowHA_ch,nHighRP_ch,nHighHA_ch, ...
    'VariableNames',{'Channel','LowCC_RP','LowCC_HA','HighCC_RP','HighCC_HA', ...
    'N_LowCC_RP','N_LowCC_HA','N_HighCC_RP','N_HighCC_HA'});
writetable(AbsoluteTopoTable,fullfile(save_dir,'Beta_FINAL_AbsoluteTopographies_SourceData_v2.csv'));

CascadeTable = table( ...
    string(chanNames(:)), ...
    CascadeLow(:,1),CascadeLow(:,2),CascadeLow(:,3),CascadeLow(:,4), ...
    CascadeHigh(:,1),CascadeHigh(:,2),CascadeHigh(:,3),CascadeHigh(:,4), ...
    'VariableNames',{'Channel','Low_200_400','Low_400_600','Low_600_800','Low_800_1000', ...
    'High_200_400','High_400_600','High_600_800','High_800_1000'});
writetable(CascadeTable,fullfile(save_dir,'Beta_FINAL_Cascade_SourceData_v2.csv'));

save(fullfile(save_dir,'Beta_FINAL_FigureSourceData_v2.mat'), ...
    'tf_file','subject_id','subject_names','idxHigh','idxLow', ...
    'times','freqs','chanNames','topoNames','topoKeep','topoLocs', ...
    'RP_roi','HA_roi','DiffROI', ...
    'meanRP','semRP','meanHA','semHA','meanLowDiff','semLowDiff','meanHighDiff','semHighDiff', ...
    'mapLowRP','mapLowHA','mapHighRP','mapHighHA','mapGrandDiff', ...
    'CascadeLow','CascadeHigh','CascadeGrand', ...
    'COMMON_CLIM','COMMON_CMAP', ...
    'display_smooth_kernel','display_interp_factor', ...
    '-v7.3');

fprintf('Figure summary\n');
fprintf('Canonical MAT: %s\n',tf_file);
fprintf('Subjects total: %d\n',N);
fprintf('High CC IDs present: %d/19\n',sum(idxHigh));
fprintf('Low CC IDs present: %d/19\n',sum(idxLow));
fprintf('Strict paired ROI9: %d/%d\n',sum(completePAIR9),N);
fprintf('Low CC paired ROI9: %d\n',sum(lowPairMask));
fprintf('High CC paired ROI9: %d\n',sum(highPairMask));
fprintf('Final topo channels: %d/%d\n',sum(topoKeep),nCh);
fprintf('Waveform display softening: triangular 3-point kernel + pchip interpolation x%d\n',display_interp_factor);
fprintf('Topoplots: style=both, contours=6, jet(256), COMMON_CLIM=[%.2f %.2f]\n',COMMON_CLIM(1),COMMON_CLIM(2));
fprintf('Output folder: %s\n',save_dir);

function [m,se,n] = local_mean_sem(X)
    n = sum(isfinite(X),2);
    m = mean(X,2,'omitnan');
    sd = std(X,0,2,'omitnan');
    se = NaN(size(m));
    good = n > 1;
    se(good) = sd(good) ./ sqrt(n(good));
end

function [m,n] = local_channel_mean(X)
    n = sum(isfinite(X),2);
    m = mean(X,2,'omitnan');
end

function [x_dense,y_dense,se_dense] = local_soft_curve(x,y,se,kernel,interp_factor)
    x = x(:);
    y = y(:);
    se = se(:);

    y_soft = local_nan_smooth(y,kernel);
    se_soft = local_nan_smooth(se,kernel);

    x_dense = linspace(x(1),x(end),numel(x)*interp_factor)';
    y_dense = interp1(x,y_soft,x_dense,'pchip');
    se_dense = interp1(x,se_soft,x_dense,'pchip');
end

function y = local_nan_smooth(x,kernel)
    x = x(:);
    good = isfinite(x);
    x0 = x;
    x0(~good) = 0;
    num = conv(x0,kernel(:),'same');
    den = conv(double(good),kernel(:),'same');
    y = num ./ den;
    y(den==0) = NaN;
end

function local_fill_sem(ax,x,m,se,col,alphaVal)
    xx = [x(:); flipud(x(:))];
    yy = [m(:)-se(:); flipud(m(:)+se(:))];
    fill(ax,xx,yy,col,'FaceAlpha',alphaVal,'EdgeColor','none','HandleVisibility','off');
end

function local_draw_topo(ax,data,chanlocs,roiIdx,climVals,cmap,showROI)
    axes(ax); %#ok<LAXES>
    cla(ax);

    if showROI && ~isempty(roiIdx)
        topoplot(data,chanlocs, ...
            'maplimits',climVals, ...
            'style','both', ...
            'numcontour',6, ...
            'electrodes','off', ...
            'shading','interp', ...
            'whitebk','on', ...
            'conv','off', ...
            'emarker2',{roiIdx,'o','k',2.4,0.7});
    else
        topoplot(data,chanlocs, ...
            'maplimits',climVals, ...
            'style','both', ...
            'numcontour',6, ...
            'electrodes','off', ...
            'shading','interp', ...
            'whitebk','on', ...
            'conv','off');
    end

    caxis(ax,climVals);
    colormap(ax,cmap);
    axis(ax,'off');

    % Use thinner topoplot lines where possible.
    lineObjs = findobj(ax,'Type','line');
    for i = 1:numel(lineObjs)
        try
            set(lineObjs(i),'LineWidth',0.6);
        catch
        end
    end
end

function cb = local_fixed_colorbar(fig,cmap,climVals,cbPosition,labelText,fontSize,labelFontSize)
    % Use a separate hidden axis for the colorbar.
    % This keeps the tiled topoplots the same size.
    axCB = axes( ...
        'Parent',fig, ...
        'Units','normalized', ...
        'Position',[0.001 0.001 0.001 0.001], ...
        'Visible','off', ...
        'HandleVisibility','off');

    colormap(axCB,cmap);
    caxis(axCB,climVals);

    cb = colorbar(axCB,'eastoutside');
    cb.Units = 'normalized';
    cb.Position = cbPosition;
    cb.Limits = climVals;
    cb.FontName = 'Arial';
    cb.FontSize = fontSize;
    cb.TickDirection = 'out';
    cb.Box = 'off';
    cb.Label.String = labelText;
    cb.Label.FontName = 'Arial';
    cb.Label.FontSize = labelFontSize;
end

function local_export_figure(fig,save_dir,baseName)
    savefig(fig,fullfile(save_dir,[baseName '.fig']));
    if exist('exportgraphics','file') == 2
        exportgraphics(fig,fullfile(save_dir,[baseName '.pdf']),'ContentType','vector');
        exportgraphics(fig,fullfile(save_dir,[baseName '.png']),'Resolution',600);
    else
        print(fig,fullfile(save_dir,[baseName '.pdf']),'-dpdf','-painters');
        print(fig,fullfile(save_dir,[baseName '.png']),'-dpng','-r600');
    end
end
