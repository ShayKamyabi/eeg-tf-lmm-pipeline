clear;
close all;
clc;
%% Paths

analysis_file = 'analysis_ready_data_extreme_groups.csv';

save_dir = 'F:\my data\Figures\Publication';
out_dir = fullfile(save_dir,'LMM1_FINAL_figures');

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

if ~isfile(analysis_file)
    [f,p] = uigetfile('*.csv','Select analysis_ready_data_extreme_groups.csv');
    if isequal(f,0)
        error('Analysis file not selected.');
    end
    analysis_file = fullfile(p,f);
end
%% Load data

data = readtable(analysis_file);

requiredVars = {'ID','CC_Group','Trial_Type','RT','Log_RT','Z_Abs_SV'};
for k = 1:numel(requiredVars)
    if ~ismember(requiredVars{k},data.Properties.VariableNames)
        error('Missing required variable: %s',requiredVars{k});
    end
end

LOW_IDS = [ ...
     1  2 11 12 13 15 27 28 38 40 ...
    45 47 48 50 52 53 60 66 67];

HIGH_IDS = [ ...
     6 17 21 22 25 29 31 32 33 36 ...
    37 39 41 46 49 54 55 58 61];

EXTREME_IDS = sort([LOW_IDS HIGH_IDS]);

ID_numeric = double(data.ID);
ids_in_file = sort(unique(ID_numeric))';

if height(data) ~= 2165
    error('Expected 2165 final LMM1 rows; found %d.',height(data));
end

if ~isequal(ids_in_file,EXTREME_IDS)
    fprintf('Expected IDs:\n'); disp(EXTREME_IDS);
    fprintf('Found IDs:\n'); disp(ids_in_file);
    error('LMM1 participant set does not match the locked final N=38 sample.');
end

if sum(strcmpi(string(data.Trial_Type),'Risky')) ~= 1482
    error('Expected 1482 risky trials.');
end

if sum(strcmpi(string(data.Trial_Type),'Ambiguous')) ~= 683
    error('Expected 683 ambiguous trials.');
end

if any(string(data.CC_Group(ismember(ID_numeric,LOW_IDS))) ~= "Low CC")
    error('At least one locked Low-CC participant is mislabeled.');
end

if any(string(data.CC_Group(ismember(ID_numeric,HIGH_IDS))) ~= "High CC")
    error('At least one locked High-CC participant is mislabeled.');
end

if any(~isfinite(data.RT)) || any(data.RT <= 0)
    error('Non-finite or non-positive RT remains in corrected LMM1 data.');
end

data.ID = categorical(data.ID);
data.CC_Group = categorical(string(data.CC_Group),{'Low CC','High CC'});
%% Split conditions

risky_data = data(strcmpi(string(data.Trial_Type),'Risky'),:);
ambig_data = data(strcmpi(string(data.Trial_Type),'Ambiguous'),:);

fprintf('Risky trials:     %d\n',height(risky_data));
fprintf('Ambiguous trials: %d\n',height(ambig_data));
fprintf('Risky subjects:   %d\n',numel(unique(risky_data.ID)));
fprintf('Ambig subjects:   %d\n',numel(unique(ambig_data.ID)));
%% Fit LMMs

formula = 'Log_RT ~ Z_Abs_SV * CC_Group + (Z_Abs_SV | ID)';

lmm_risky = fitlme(risky_data,formula,'FitMethod','REML');
lmm_ambig = fitlme(ambig_data,formula,'FitMethod','REML');

anova_risky = anova(lmm_risky);
anova_ambig = anova(lmm_ambig);

disp(anova_risky);

disp(anova_ambig);
%% Effect estimates

statsRisky = extract_lmm1_stats(lmm_risky,anova_risky);
statsAmbig = extract_lmm1_stats(lmm_ambig,anova_ambig);

fprintf('Risky:     b = %.6f, 95%% CI [%.6f, %.6f], p = %.8f\n', ...
    statsRisky.Interaction_b,statsRisky.Interaction_CI(1), ...
    statsRisky.Interaction_CI(2),statsRisky.Interaction_p);
fprintf('Ambiguous: b = %.6f, 95%% CI [%.6f, %.6f], p = %.8f\n', ...
    statsAmbig.Interaction_b,statsAmbig.Interaction_CI(1), ...
    statsAmbig.Interaction_CI(2),statsAmbig.Interaction_p);
%% Plot settings

COL_LOW  = [65 65 65]/255;
COL_HIGH = [0 114 178]/255;

COL_LOW_FILL  = [0.72 0.72 0.72];
COL_HIGH_FILL = [0.45 0.68 0.82];

COL_RISKY = [76 114 176]/255;
COL_AMBIG = [202 83 66]/255;

COL_TEXT      = [35 35 35]/255;
COL_GREY      = [115 115 115]/255;
COL_LIGHTGREY = [0.94 0.94 0.94];
COL_ZERO      = [0.68 0.68 0.68];
COL_GRID      = [0.90 0.90 0.90];

FONT_NAME = 'Arial';
nBins = 5;
%% Common support

[risky_zmin,risky_zmax] = common_support(risky_data);
[ambig_zmin,ambig_zmax] = common_support(ambig_data);

fprintf('Risky:     %.4f to %.4f\n',risky_zmin,risky_zmax);
fprintf('Ambiguous: %.4f to %.4f\n',ambig_zmin,ambig_zmax);
%% Model predictions

predRisky = make_predictions(lmm_risky,risky_data,risky_zmin,risky_zmax);
predAmbig = make_predictions(lmm_ambig,ambig_data,ambig_zmin,ambig_zmax);
%% Shared x-axis range

sharedConflictMin = max([min(predRisky.Display_Conflict), min(predAmbig.Display_Conflict)]);
sharedConflictMax = min([max(predRisky.Display_Conflict), max(predAmbig.Display_Conflict)]);

if sharedConflictMin >= sharedConflictMax
    error('Risky and Ambiguous conditions have no common display support.');
end

fprintf('%.4f to %.4f\n',sharedConflictMin,sharedConflictMax);
fprintf('RIGHTWARD = GREATER CONFLICT\n');
%% Observed bins

[summaryRisky,subjectBinsRisky] = make_observed_bins(risky_data,nBins,risky_zmin,risky_zmax);
[summaryAmbig,subjectBinsAmbig] = make_observed_bins(ambig_data,nBins,ambig_zmin,ambig_zmax);
%% Simple slopes

slopesRisky = get_display_conflict_slopes(lmm_risky);
slopesAmbig = get_display_conflict_slopes(lmm_ambig);
%% Shared y-axis range

allY = [ ...
    predRisky.CI_Lower; predRisky.CI_Upper; ...
    predAmbig.CI_Lower; predAmbig.CI_Upper; ...
    summaryRisky.Mean_Log_RT-summaryRisky.SEM_Log_RT; ...
    summaryRisky.Mean_Log_RT+summaryRisky.SEM_Log_RT; ...
    summaryAmbig.Mean_Log_RT-summaryAmbig.SEM_Log_RT; ...
    summaryAmbig.Mean_Log_RT+summaryAmbig.SEM_Log_RT];

allY = allY(isfinite(allY));
yMin = min(allY);
yMax = max(allY);
yRange = yMax-yMin;
if yRange <= 0, yRange = 1; end
sharedYLim = [yMin-0.09*yRange yMax+0.09*yRange];
%% Interaction contrasts

displayRisky_b  = -statsRisky.Interaction_b;
displayRisky_CI = [-statsRisky.Interaction_CI(2), -statsRisky.Interaction_CI(1)];

displayAmbig_b  = -statsAmbig.Interaction_b;
displayAmbig_CI = [-statsAmbig.Interaction_CI(2), -statsAmbig.Interaction_CI(1)];

fprintf('Displayed conflict = -Z(|SV|)\n');
fprintf('Risky:     %.4f [%.4f, %.4f], p = %.6f\n', ...
    displayRisky_b,displayRisky_CI(1),displayRisky_CI(2),statsRisky.Interaction_p);
fprintf('Ambiguous: %.4f [%.4f, %.4f], p = %.6f\n', ...
    displayAmbig_b,displayAmbig_CI(1),displayAmbig_CI(2),statsAmbig.Interaction_p);
%% Effect arrays

slopeEffect = [ ...
    slopesRisky.Low.Estimate; ...
    slopesRisky.High.Estimate; ...
    slopesAmbig.Low.Estimate; ...
    slopesAmbig.High.Estimate];

slopeLower = [ ...
    slopesRisky.Low.CI_Lower; ...
    slopesRisky.High.CI_Lower; ...
    slopesAmbig.Low.CI_Lower; ...
    slopesAmbig.High.CI_Lower];

slopeUpper = [ ...
    slopesRisky.Low.CI_Upper; ...
    slopesRisky.High.CI_Upper; ...
    slopesAmbig.Low.CI_Upper; ...
    slopesAmbig.High.CI_Upper];

contrastEffect = [displayRisky_b; displayAmbig_b];
contrastLower  = [displayRisky_CI(1); displayAmbig_CI(1)];
contrastUpper  = [displayRisky_CI(2); displayAmbig_CI(2)];
contrastP      = [statsRisky.Interaction_p; statsAmbig.Interaction_p];
%% Figure A

figA = figure( ...
    'Name','LMM1 A - Conflict-index definition', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 15.8 8.8], ...
    'Renderer','painters');

axA = axes( ...
    'Parent',figA, ...
    'Units','normalized', ...
    'Position',[0.05 0.08 0.90 0.80]);

hold(axA,'on');
axis(axA,'off');
xlim(axA,[0 1]);
ylim(axA,[0 1]);

draw_info_box_q1(axA,[0.04 0.38 0.25 0.31], ...
    [0.97 0.97 0.97],[0.55 0.55 0.55], ...
    '|SV|', ...
    {'smaller |SV|';'closer to indifference';'greater conflict'}, ...
    COL_TEXT);

quiver(axA,0.305,0.535,0.09,0,0, ...
    'Color',COL_GREY, ...
    'LineWidth',1.1, ...
    'MaxHeadSize',0.22);

draw_info_box_q1(axA,[0.39 0.38 0.25 0.31], ...
    [0.965 0.98 1.00],COL_HIGH, ...
    'Fitted predictor', ...
    {'z(|SV|)';'0 = mean |SV|';'not peak conflict'}, ...
    COL_TEXT);

quiver(axA,0.655,0.535,0.09,0,0, ...
    'Color',COL_GREY, ...
    'LineWidth',1.1, ...
    'MaxHeadSize',0.22);

draw_info_box_q1(axA,[0.74 0.38 0.22 0.31], ...
    [0.965 0.98 1.00],COL_HIGH, ...
    'Display only', ...
    {'-z(|SV|)';'rightward';'= greater conflict'}, ...
    COL_TEXT);

text(axA,0.50,0.22, ...
    'All LMMs are fitted on +z(|SV|); sign reversal is used only for figure display.', ...
    'HorizontalAlignment','center', ...
    'FontName',FONT_NAME, ...
    'FontSize',8.2, ...
    'Color',COL_GREY);

title(axA, ...
    'A   Conflict-index definition', ...
    'FontName',FONT_NAME, ...
    'FontSize',11.5, ...
    'FontWeight','normal', ...
    'Color',COL_TEXT);

export_one_figure( ...
    figA,out_dir, ...
    'LMM1_A_ConflictIndexDefinition_Q1', ...
    15.8,8.8);
%% Figure B

figB = figure( ...
    'Name','LMM1 B - Risky decisions', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 12.8 10.2], ...
    'Renderer','painters');

axB = axes( ...
    'Parent',figB, ...
    'Units','normalized', ...
    'Position',[0.16 0.20 0.77 0.68]);

plot_condition_panel_q1( ...
    axB,predRisky,summaryRisky, ...
    COL_LOW,COL_HIGH,COL_LOW_FILL,COL_HIGH_FILL,COL_GRID);

xlim(axB,[sharedConflictMin sharedConflictMax]);
ylim(axB,sharedYLim);

xlabel(axB, ...
    'Model-derived conflict index,  -z(|SV|)', ...
    'FontName',FONT_NAME, ...
    'FontSize',9.4);

ylabel(axB, ...
    'Log reaction time', ...
    'FontName',FONT_NAME, ...
    'FontSize',9.4);

title(axB, ...
    'B   Risky decisions', ...
    'FontName',FONT_NAME, ...
    'FontSize',11.5, ...
    'FontWeight','normal');

annotation(figB,'textbox',[0.16 0.075 0.28 0.035], ...
    'String','lower conflict', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','left', ...
    'FontName',FONT_NAME, ...
    'FontSize',7.6, ...
    'Color',COL_GREY);

annotation(figB,'textbox',[0.65 0.075 0.28 0.035], ...
    'String','greater conflict', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','right', ...
    'FontName',FONT_NAME, ...
    'FontSize',7.6, ...
    'FontWeight','bold', ...
    'Color',COL_GREY);

export_one_figure( ...
    figB,out_dir, ...
    'LMM1_B_RiskyDecisions_Q1', ...
    12.8,10.2);
%% Figure C

figC = figure( ...
    'Name','LMM1 C - Ambiguous decisions', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[4 4 12.8 10.2], ...
    'Renderer','painters');

axC = axes( ...
    'Parent',figC, ...
    'Units','normalized', ...
    'Position',[0.16 0.20 0.77 0.68]);

plot_condition_panel_q1( ...
    axC,predAmbig,summaryAmbig, ...
    COL_LOW,COL_HIGH,COL_LOW_FILL,COL_HIGH_FILL,COL_GRID);

xlim(axC,[sharedConflictMin sharedConflictMax]);
ylim(axC,sharedYLim);

xlabel(axC, ...
    'Model-derived conflict index,  -z(|SV|)', ...
    'FontName',FONT_NAME, ...
    'FontSize',9.4);

ylabel(axC, ...
    'Log reaction time', ...
    'FontName',FONT_NAME, ...
    'FontSize',9.4);

title(axC, ...
    'C   Ambiguous decisions', ...
    'FontName',FONT_NAME, ...
    'FontSize',11.5, ...
    'FontWeight','normal');

annotation(figC,'textbox',[0.16 0.075 0.28 0.035], ...
    'String','lower conflict', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','left', ...
    'FontName',FONT_NAME, ...
    'FontSize',7.6, ...
    'Color',COL_GREY);

annotation(figC,'textbox',[0.65 0.075 0.28 0.035], ...
    'String','greater conflict', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','right', ...
    'FontName',FONT_NAME, ...
    'FontSize',7.6, ...
    'FontWeight','bold', ...
    'Color',COL_GREY);

export_one_figure( ...
    figC,out_dir, ...
    'LMM1_C_AmbiguousDecisions_Q1', ...
    12.8,10.2);
%% Figure D

figD = figure( ...
    'Name','LMM1 D - Model-estimated conflict sensitivity', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[5 5 15.8 10.4], ...
    'Renderer','painters');

axD = axes( ...
    'Parent',figD, ...
    'Units','normalized', ...
    'Position',[0.16 0.18 0.78 0.68]);

hold(axD,'on');

xRisk = [1-0.11 1+0.11];
xAmb  = [2-0.11 2+0.11];
xPts = [xRisk xAmb];

cols = [COL_LOW;COL_HIGH;COL_LOW;COL_HIGH];

for k = 1:4
    errLow = slopeEffect(k)-slopeLower(k);
    errUp  = slopeUpper(k)-slopeEffect(k);

    errorbar(axD, ...
        xPts(k),slopeEffect(k),errLow,errUp,'o', ...
        'Color',cols(k,:), ...
        'MarkerFaceColor',cols(k,:), ...
        'MarkerEdgeColor','w', ...
        'MarkerSize',7.4, ...
        'LineWidth',1.45, ...
        'CapSize',5, ...
        'HandleVisibility','off');
end

yline(axD,0,'--', ...
    'Color',COL_ZERO, ...
    'LineWidth',0.9, ...
    'HandleVisibility','off');

hLow = plot(axD,nan,nan,'o', ...
    'MarkerFaceColor',COL_LOW, ...
    'MarkerEdgeColor','w', ...
    'MarkerSize',7.4, ...
    'Color',COL_LOW, ...
    'LineWidth',1.2);

hHigh = plot(axD,nan,nan,'o', ...
    'MarkerFaceColor',COL_HIGH, ...
    'MarkerEdgeColor','w', ...
    'MarkerSize',7.4, ...
    'Color',COL_HIGH, ...
    'LineWidth',1.2);

set(axD, ...
    'XTick',[1 2], ...
    'XTickLabel',{'Risky','Ambiguous'});

xlim(axD,[0.55 2.45]);

ylabel(axD, ...
    'Change in log RT per +1 SD displayed conflict', ...
    'FontName',FONT_NAME, ...
    'FontSize',9.1);

title(axD, ...
    'D   Model-estimated conflict sensitivity', ...
    'FontName',FONT_NAME, ...
    'FontSize',11.5, ...
    'FontWeight','normal');

legend(axD,[hLow hHigh],{'Low CC','High CC'}, ...
    'Location','northwest', ...
    'Box','off', ...
    'FontName',FONT_NAME, ...
    'FontSize',8.0);

set_pub_axes(axD);
axD.YGrid = 'on';
axD.GridColor = COL_GRID;
axD.GridAlpha = 0.24;

local_pad_ylim(axD,0.13);

export_one_figure( ...
    figD,out_dir, ...
    'LMM1_D_ConflictSensitivity_Q1', ...
    15.8,10.4);
%% Figure E

figE = figure( ...
    'Name','LMM1 E - CC-group interaction contrast', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[6 6 13.4 10.2], ...
    'Renderer','painters');

axE = axes( ...
    'Parent',figE, ...
    'Units','normalized', ...
    'Position',[0.18 0.18 0.75 0.68]);

hold(axE,'on');

for k = 1:2
    errLow = contrastEffect(k)-contrastLower(k);
    errUp  = contrastUpper(k)-contrastEffect(k);

    if k == 1
        thisCol = COL_RISKY;
    else
        thisCol = COL_AMBIG;
    end

    errorbar(axE, ...
        k,contrastEffect(k),errLow,errUp,'o', ...
        'Color',thisCol, ...
        'MarkerFaceColor',thisCol, ...
        'MarkerEdgeColor','w', ...
        'MarkerSize',7.8, ...
        'LineWidth',1.55, ...
        'CapSize',5, ...
        'HandleVisibility','off');
end

yline(axE,0,'--', ...
    'Color',COL_ZERO, ...
    'LineWidth',0.9, ...
    'HandleVisibility','off');

set(axE, ...
    'XTick',[1 2], ...
    'XTickLabel',{'Risky','Ambiguous'});

xlim(axE,[0.55 2.45]);

ylabel(axE, ...
    'High - Low slope contrast', ...
    'FontName',FONT_NAME, ...
    'FontSize',9.1);

title(axE, ...
    'E   CC-group difference in conflict sensitivity', ...
    'FontName',FONT_NAME, ...
    'FontSize',11.5, ...
    'FontWeight','normal');

set_pub_axes(axE);
axE.YGrid = 'on';
axE.GridColor = COL_GRID;
axE.GridAlpha = 0.24;

local_pad_ylim(axE,0.25);

yl = ylim(axE);
yr = diff(yl);
if yr <= 0, yr = 1; end

for k = 1:2
    yText = contrastUpper(k) + 0.055*yr;
    text(axE,k,yText,format_p(contrastP(k)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontName',FONT_NAME, ...
        'FontSize',8.0, ...
        'Color',COL_TEXT);
end

export_one_figure( ...
    figE,out_dir, ...
    'LMM1_E_GroupInteractionContrast_Q1', ...
    13.4,10.2);

fprintf('Output folder:\n%s\n',out_dir);
%% Export data

writetable(predRisky,fullfile(out_dir,'LMM1_Risky_Predictions.csv'));
writetable(predAmbig,fullfile(out_dir,'LMM1_Ambiguous_Predictions.csv'));
writetable(summaryRisky,fullfile(out_dir,'LMM1_Risky_ObservedBins.csv'));
writetable(summaryAmbig,fullfile(out_dir,'LMM1_Ambiguous_ObservedBins.csv'));
writetable(subjectBinsRisky,fullfile(out_dir,'LMM1_Risky_SubjectBins.csv'));
writetable(subjectBinsAmbig,fullfile(out_dir,'LMM1_Ambiguous_SubjectBins.csv'));

SlopeTable = table( ...
    ["Risky";"Risky";"Ambiguous";"Ambiguous"], ...
    ["Low CC";"High CC";"Low CC";"High CC"], ...
    slopeEffect,slopeLower,slopeUpper, ...
    'VariableNames',{'Condition','CC_Group','DisplayedConflictSlope','CI_Lower','CI_Upper'});
writetable(SlopeTable,fullfile(out_dir,'LMM1_DisplayConflictSlopes.csv'));

InteractionTable = table( ...
    ["Risky";"Ambiguous"], ...
    contrastEffect,contrastLower,contrastUpper,contrastP, ...
    'VariableNames',{'Condition','HighMinusLow_DisplaySlope','Lower95','Upper95','pValue'});
writetable(InteractionTable,fullfile(out_dir,'LMM1_DisplayInteractionContrasts.csv'));

save(fullfile(out_dir,'LMM1_FINAL_Q1_STYLE_SEPARATE_RESULTS.mat'), ...
    'lmm_risky','lmm_ambig','anova_risky','anova_ambig', ...
    'statsRisky','statsAmbig','slopesRisky','slopesAmbig', ...
    'predRisky','predAmbig','summaryRisky','summaryAmbig', ...
    'subjectBinsRisky','subjectBinsAmbig','SlopeTable','InteractionTable');
%% Common support

function [zMin,zMax] = common_support(D)
    low = D.Z_Abs_SV(string(D.CC_Group)=="Low CC");
    high = D.Z_Abs_SV(string(D.CC_Group)=="High CC");

    lowRange = [prctile(low,2.5) prctile(low,97.5)];
    highRange = [prctile(high,2.5) prctile(high,97.5)];

    zMin = max(lowRange(1),highRange(1));
    zMax = min(lowRange(2),highRange(2));

    if zMin >= zMax
        error('No usable overlapping Z|SV| support between CC groups.');
    end
end
%% Predictions

function T = make_predictions(model,D,zMin,zMax)
    zGrid = linspace(zMin,zMax,220)';
    groupNames = ["Low CC","High CC"];
    T = table();

    for g = 1:2
        n = numel(zGrid);
        groupVar = categorical(repmat(groupNames(g),n,1),{'Low CC','High CC'});
        idVar = repmat(D.ID(1),n,1);
        Tpred = table(groupVar,zGrid,idVar,'VariableNames',{'CC_Group','Z_Abs_SV','ID'});

        [yPred,yCI] = predict(model,Tpred,'Conditional',false);

        Tout = table();
        Tout.CC_Group = repmat(groupNames(g),n,1);
        Tout.Display_Conflict = -zGrid;
        Tout.Predicted_Log_RT = yPred;
        Tout.CI_Lower = yCI(:,1);
        Tout.CI_Upper = yCI(:,2);

        T = [T;Tout]; %#ok<AGROW>
    end
end
%% Observed bins

function [Summary,SubjectBins] = make_observed_bins(D,nBins,zMin,zMax)
    D = D(D.Z_Abs_SV >= zMin & D.Z_Abs_SV <= zMax,:);
    pooledZ = D.Z_Abs_SV;

    edges = quantile(pooledZ,linspace(0,1,nBins+1));
    edges = unique(edges,'stable');

    if numel(edges) < 3
        error('Insufficient unique Z|SV| values to construct bins.');
    end

    edges(1) = min(pooledZ)-1e-10;
    edges(end) = max(pooledZ)+1e-10;

    IDs = unique(D.ID);

    SID = strings(0,1);
    GROUP = strings(0,1);
    BIN = zeros(0,1);
    CONFLICT = zeros(0,1);
    LOGRT = zeros(0,1);
    NTRIALS = zeros(0,1);

    for s = 1:numel(IDs)
        Ds = D(D.ID==IDs(s),:);
        bins = discretize(Ds.Z_Abs_SV,edges);
        validBins = unique(bins(~isnan(bins)));

        for bb = 1:numel(validBins)
            b = validBins(bb);
            idx = bins==b;
            if ~any(idx), continue; end

            SID(end+1,1) = string(IDs(s)); %#ok<AGROW>
            GROUP(end+1,1) = string(Ds.CC_Group(1)); %#ok<AGROW>
            BIN(end+1,1) = b; %#ok<AGROW>
            CONFLICT(end+1,1) = -mean(Ds.Z_Abs_SV(idx),'omitnan'); %#ok<AGROW>
            LOGRT(end+1,1) = mean(Ds.Log_RT(idx),'omitnan'); %#ok<AGROW>
            NTRIALS(end+1,1) = sum(idx); %#ok<AGROW>
        end
    end

    SubjectBins = table(SID,GROUP,BIN,CONFLICT,LOGRT,NTRIALS, ...
        'VariableNames',{'Subject','CC_Group','Bin','Display_Conflict','Log_RT','N_Trials'});

    groupNames = ["Low CC","High CC"];

    OUT_GROUP = strings(0,1);
    OUT_BIN = zeros(0,1);
    MEAN_X = zeros(0,1);
    MEAN_Y = zeros(0,1);
    SEM_Y = zeros(0,1);
    NSUB = zeros(0,1);

    allBins = unique(SubjectBins.Bin);

    for g = 1:numel(groupNames)
        for bb = 1:numel(allBins)
            b = allBins(bb);
            idx = SubjectBins.CC_Group==groupNames(g) & SubjectBins.Bin==b;
            if ~any(idx), continue; end

            OUT_GROUP(end+1,1) = groupNames(g); %#ok<AGROW>
            OUT_BIN(end+1,1) = b; %#ok<AGROW>
            MEAN_X(end+1,1) = mean(SubjectBins.Display_Conflict(idx),'omitnan'); %#ok<AGROW>
            MEAN_Y(end+1,1) = mean(SubjectBins.Log_RT(idx),'omitnan'); %#ok<AGROW>

            n = sum(idx);
            NSUB(end+1,1) = n; %#ok<AGROW>

            if n > 1
                SEM_Y(end+1,1) = std(SubjectBins.Log_RT(idx),'omitnan')/sqrt(n); %#ok<AGROW>
            else
                SEM_Y(end+1,1) = NaN; %#ok<AGROW>
            end
        end
    end

    Summary = table(OUT_GROUP,OUT_BIN,MEAN_X,MEAN_Y,SEM_Y,NSUB, ...
        'VariableNames',{'CC_Group','Bin','Mean_Display_Conflict','Mean_Log_RT','SEM_Log_RT','N_Subjects'});
end
%% Simple slopes

function S = get_display_conflict_slopes(model)
    [idxZ,idxInt] = get_lmm1_coef_indices(model);

    C = model.Coefficients;
    b = C.Estimate;
    V = model.CoefficientCovariance;
    df = model.DFE;
    tCrit = tinv(0.975,df);

    lowEstimate = -b(idxZ);
    lowSE = sqrt(V(idxZ,idxZ));
    lowCI = [lowEstimate-tCrit*lowSE, lowEstimate+tCrit*lowSE];

    highEstimate = -(b(idxZ)+b(idxInt));
    highVar = V(idxZ,idxZ)+V(idxInt,idxInt)+2*V(idxZ,idxInt);
    highSE = sqrt(max(highVar,0));
    highCI = [highEstimate-tCrit*highSE, highEstimate+tCrit*highSE];

    S.Low.Estimate = lowEstimate;
    S.Low.SE = lowSE;
    S.Low.CI_Lower = lowCI(1);
    S.Low.CI_Upper = lowCI(2);

    S.High.Estimate = highEstimate;
    S.High.SE = highSE;
    S.High.CI_Lower = highCI(1);
    S.High.CI_Upper = highCI(2);
end
%% Prediction plot

function plot_condition_panel_q1(ax,Pred,Summary,colLow,colHigh,fillLow,fillHigh,colGrid)
    hold(ax,'on');

    groupNames = ["Low CC","High CC"];
    lineCols = [colLow;colHigh];
    fillCols = [fillLow;fillHigh];
    fillAlpha = [0.38 0.32];
    hLine = gobjects(2,1);

    for g = 1:2
        P = Pred(Pred.CC_Group==groupNames(g),:);
        [~,ord] = sort(P.Display_Conflict);
        P = P(ord,:);

        fill(ax,[P.Display_Conflict;flipud(P.Display_Conflict)], ...
            [P.CI_Lower;flipud(P.CI_Upper)], ...
            fillCols(g,:), ...
            'FaceAlpha',fillAlpha(g), ...
            'EdgeColor','none', ...
            'HandleVisibility','off');
    end

    for g = 1:2
        P = Pred(Pred.CC_Group==groupNames(g),:);
        [~,ord] = sort(P.Display_Conflict);
        P = P(ord,:);

        hLine(g) = plot(ax,P.Display_Conflict,P.Predicted_Log_RT,'-', ...
            'Color',lineCols(g,:),'LineWidth',2.45);
    end

    for g = 1:2
        O = Summary(Summary.CC_Group==groupNames(g),:);
        if isempty(O), continue; end

        errorbar(ax,O.Mean_Display_Conflict,O.Mean_Log_RT,O.SEM_Log_RT,'o', ...
            'Color',lineCols(g,:), ...
            'MarkerFaceColor','w', ...
            'MarkerEdgeColor',lineCols(g,:), ...
            'MarkerSize',5.2, ...
            'LineWidth',1.00, ...
            'CapSize',3.2, ...
            'HandleVisibility','off');
    end

    legend(ax,hLine,{'Low CC','High CC'}, ...
        'Location','northwest','Box','off','FontName','Arial','FontSize',7.7);

    set_pub_axes(ax);
    ax.YGrid = 'on';
    ax.XGrid = 'off';
    ax.GridColor = colGrid;
    ax.GridAlpha = 0.24;
end
%% LMM1 statistics

function S = extract_lmm1_stats(model,A)
    [b,~,coefStats] = fixedEffects(model);
    ci = coefCI(model);

    if numel(b) < 4
        error('Expected at least 4 fixed-effect coefficients in LMM1.');
    end

    S = struct();
    S.Interaction_b = b(4);
    S.Interaction_CI = [ci(4,1) ci(4,2)];

    try
        S.Interaction_SE = coefStats.SE(4);
    catch
        V = model.CoefficientCovariance;
        S.Interaction_SE = sqrt(V(4,4));
    end

    try
        S.Interaction_p = coefStats.pValue(4);
    catch
        S.Interaction_p = NaN;
    end

    terms = get_anova_terms(A);
    idx = [];

    for ii = 1:numel(terms)
        term_i = char(terms(ii));
        if ~isempty(strfind(term_i,'CC_Group')) && ~isempty(strfind(term_i,'Z_Abs_SV')) %#ok<STREMP>
            idx(end+1) = ii; %#ok<AGROW>
        end
    end

    if numel(idx) ~= 1
        fprintf('\nANOVA terms found:\n');
        disp(terms);
        error('Could not uniquely identify CC_Group x Z_Abs_SV interaction.');
    end

    try
        S.Interaction_F = double(A.FStat(idx));
        S.Interaction_DF1 = double(A.DF1(idx));
        S.Interaction_DF2 = double(A.DF2(idx));
        S.Interaction_p = double(A.pValue(idx));
    catch
        S.Interaction_F = NaN;
        S.Interaction_DF1 = NaN;
        S.Interaction_DF2 = NaN;
    end
end

function [idxZ,idxInt] = get_lmm1_coef_indices(model)
    names = string(model.Coefficients.Name);
    idxZ = find(names=="Z_Abs_SV");

    if numel(idxZ) ~= 1
        error('Could not identify Z_Abs_SV coefficient.');
    end

    normalized = erase(names,[" ","(",")","{","}","'"]);
    maskInt = contains(normalized,'CC_Group') & ...
        contains(normalized,'Z_Abs_SV') & ...
        count(normalized,":")==1;

    idxInt = find(maskInt);

    if numel(idxInt) ~= 1
        fprintf('\nCoefficient names found:\n');
        disp(names);
        error('Could not uniquely identify CC_Group x Z_Abs_SV coefficient.');
    end
end

function terms = get_anova_terms(A)
    try
        terms = string(A.Term);
        terms = strip(terms(:));
        if ~isempty(terms), return; end
    catch
    end

    if istable(A) && ~isempty(A.Properties.RowNames)
        terms = string(A.Properties.RowNames);
        terms = strip(terms(:));
        return;
    end

    error('Could not determine ANOVA term names.');
end
%% Schematic box

function draw_info_box_q1(ax,pos,faceCol,edgeCol,header,body,textCol)
    rectangle(ax,'Position',pos,'Curvature',0.035, ...
        'FaceColor',faceCol,'EdgeColor',edgeCol,'LineWidth',1.0);

    cx = pos(1)+pos(3)/2;

    text(ax,cx,pos(2)+0.185,header, ...
        'HorizontalAlignment','center', ...
        'FontName','Arial','FontSize',8.8,'FontWeight','bold','Color',edgeCol);

    text(ax,cx,pos(2)+0.090,body, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Arial','FontSize',7.2,'Color',textCol);
end
%% Conflict labels

function local_conflict_direction_labels(ax,colGrey,fontName)
    text(ax,0.01,-0.18,'lower conflict', ...
        'Units','normalized','HorizontalAlignment','left', ...
        'FontName',fontName,'FontSize',7.0,'Color',colGrey,'Clipping','off');

    text(ax,0.99,-0.18,'greater conflict', ...
        'Units','normalized','HorizontalAlignment','right', ...
        'FontName',fontName,'FontSize',7.0,'FontWeight','bold', ...
        'Color',colGrey,'Clipping','off');
end
%% Axes

function set_pub_axes(ax)
    set(ax,'FontName','Arial','FontSize',8.2,'LineWidth',0.85, ...
        'TickDir','out','TickLength',[0.018 0.018],'Box','off','Layer','top');
end

function local_pad_ylim(ax,fraction)
    yl = ylim(ax);
    r = diff(yl);
    if r <= 0, r = 1; end
    ylim(ax,[yl(1)-fraction*r yl(2)+fraction*r]);
end
%% Figure export

function export_one_figure(fig,out_dir,baseName,w_cm,h_cm)
    drawnow;
    pause(0.20);

    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    testFile = fullfile(out_dir,'__matlab_write_test__.tmp');
    fid = fopen(testFile,'w');
    if fid < 0
        error('Output folder is not writable: %s',out_dir);
    end
    fclose(fid);
    delete(testFile);

    fig_file = fullfile(out_dir,[baseName '.fig']);
    pdf_file = fullfile(out_dir,[baseName '.pdf']);
    eps_file = fullfile(out_dir,[baseName '.eps']);
    png_file = fullfile(out_dir,[baseName '.png']);
    tif_file = fullfile(out_dir,[baseName '.tif']);
    svg_file = fullfile(out_dir,[baseName '.svg']);

    safe_delete(pdf_file);
    safe_delete(eps_file);
    safe_delete(png_file);
    safe_delete(tif_file);
    safe_delete(svg_file);

    savefig(fig,fig_file);

    set(fig,'PaperUnits','centimeters');
    set(fig,'PaperSize',[w_cm h_cm]);
    set(fig,'PaperPosition',[0 0 w_cm h_cm]);
    set(fig,'PaperPositionMode','manual');
    set(fig,'InvertHardcopy','off');

    print(fig,pdf_file,'-dpdf','-painters');
    print(fig,eps_file,'-depsc','-painters');

    temp_png = fullfile(tempdir,[baseName '__temp.png']);
    temp_tif = fullfile(tempdir,[baseName '__temp.tif']);
    safe_delete(temp_png);
    safe_delete(temp_tif);

    try
        exportgraphics(fig,temp_png,'Resolution',600,'BackgroundColor','white');
        [ok,msg] = copyfile(temp_png,png_file,'f');
        if ~ok, error('Could not copy PNG to destination: %s',msg); end
    catch MEpng
        warning('PNG export failed: %s',MEpng.message);
    end

    try
        exportgraphics(fig,temp_tif,'Resolution',600,'BackgroundColor','white');
        [ok,msg] = copyfile(temp_tif,tif_file,'f');
        if ~ok, error('Could not copy TIFF to destination: %s',msg); end
    catch MEtif
        warning('TIFF export failed: %s',MEtif.message);
    end

    safe_delete(temp_png);
    safe_delete(temp_tif);

    try
        exportgraphics(fig,svg_file,'ContentType','vector','BackgroundColor','white');
    catch MEsvg
        warning('SVG export skipped/failed: %s',MEsvg.message);
    end

    fprintf('\nExported: %s\n',baseName);
    fprintf('  FIG : %s\n',fig_file);
    fprintf('  PDF : %s\n',pdf_file);
    fprintf('  EPS : %s\n',eps_file);
    fprintf('  PNG : %s\n',png_file);
    fprintf('  TIF : %s\n',tif_file);
end

function safe_delete(filePath)
    if isfile(filePath)
        try
            delete(filePath);
        catch
            warning('Could not delete existing file before export: %s',filePath);
        end
    end
end
%% P-value format

function txt = format_p(p)
    if p < 0.001
        txt = 'p < .001';
    else
        txt = sprintf('p = %.3f',p);
    end
end
