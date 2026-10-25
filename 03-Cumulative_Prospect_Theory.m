clear;
close all;
clc;

save_dir = 'F:\my data\Figures\Publication';

params_file = fullfile(save_dir,'CPT_SubjectParameters_N57.csv');
trial_file  = fullfile(save_dir,'CPT_TrialLevel_N57.csv');

out_dir = fullfile(save_dir,'CPT_SEPARATE_REFINED_FINAL');

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

if ~exist(params_file,'file')
    error('Cannot find: %s', params_file);
end

if ~exist(trial_file,'file')
    error('Cannot find: %s', trial_file);
end

P = readtable(params_file);
T = readtable(trial_file);

if ismember('ExitFlag', P.Properties.VariableNames)
    P = P(P.ExitFlag > 0, :);
end

fprintf('CPT figures\n');
fprintf('Participants in parameter table: %d\n', height(P));
fprintf('Unique participants in trial table: %d\n', numel(unique(string(T.ID))));

med_alpha  = median(P.alpha,  'omitnan');
med_lambda = median(P.lambda, 'omitnan');
med_gamma  = median(P.gamma,  'omitnan');
med_tau    = median(P.tau,    'omitnan'); %#ok<NASGU>

accRisky = P.Accuracy_Risky;
accAmbig = P.Accuracy_Ambiguous;
r2Risky  = P.McFadden_R2_Risky;

pGrid = linspace(0.001,0.999,400)';
wGrid = CPT_weight(pGrid, med_gamma);

outcomeGrid = linspace(-20,20,500)';
utilityGrid = zeros(size(outcomeGrid));

idxGain = outcomeGrid >= 0;
idxLoss = outcomeGrid < 0;

utilityGrid(idxGain) = outcomeGrid(idxGain).^med_alpha;
utilityGrid(idxLoss) = -med_lambda .* abs(outcomeGrid(idxLoss)).^med_alpha;

COL_RISKY = [0 114 178] / 255;
COL_AMBIG = [213 94 0] / 255;
COL_CPT   = [0 158 115] / 255;

COL_DARK    = [45 45 45] / 255;
COL_GREY    = [120 120 120] / 255;
COL_MIDGREY = [175 175 175] / 255;

COL_RISKY_LIGHT = 0.65 + 0.35*COL_RISKY;
COL_AMBIG_LIGHT = 0.65 + 0.35*COL_AMBIG;
COL_CPT_LIGHT   = 0.70 + 0.30*COL_CPT;

set(0,'DefaultAxesFontName','Arial');

figA = figure( ...
    'Name','CPT A - Computational architecture', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 16.0 11.0], ...
    'Renderer','painters');

axA = axes('Parent',figA,'Units','normalized','Position',[0.05 0.08 0.90 0.84]);
hold(axA,'on');
axis(axA,'off');
xlim(axA,[0 1]);
ylim(axA,[0 1]);

rectangle(axA,'Position',[0.04 0.64 0.24 0.20],'Curvature',0.06, ...
    'FaceColor',[0.96 0.98 1.00],'EdgeColor',COL_RISKY,'LineWidth',1.3);
text(axA,0.16,0.765,'Risky','HorizontalAlignment','center', ...
    'FontSize',10,'FontWeight','bold','Color',COL_RISKY);
text(axA,0.16,0.695,{'Gain / loss';'Known probabilities'}, ...
    'HorizontalAlignment','center','FontSize',8.2,'Color',COL_DARK);

rectangle(axA,'Position',[0.04 0.24 0.24 0.24],'Curvature',0.06, ...
    'FaceColor',[1.00 0.97 0.95],'EdgeColor',COL_AMBIG,'LineWidth',1.3);
text(axA,0.16,0.405,'Ambiguous','HorizontalAlignment','center', ...
    'FontSize',10,'FontWeight','bold','Color',COL_AMBIG);
text(axA,0.16,0.325,{'Missing information';'Model-coded for valuation'}, ...
    'HorizontalAlignment','center','FontSize',8.0,'Color',COL_DARK);
text(axA,0.16,0.265,'p = 0.5 for hidden probability', ...
    'HorizontalAlignment','center','FontSize',7.1,'Color',COL_GREY);

rectangle(axA,'Position',[0.38 0.40 0.24 0.28],'Curvature',0.06, ...
    'FaceColor',[0.95 0.99 0.98],'EdgeColor',COL_CPT,'LineWidth',1.4);
text(axA,0.50,0.605,'CPT valuation','HorizontalAlignment','center', ...
    'FontSize',10,'FontWeight','bold','Color',COL_CPT);
text(axA,0.50,0.525,{'Outcome utility: \alpha, \lambda';'Probability weighting: \gamma'}, ...
    'HorizontalAlignment','center','FontSize',8.2,'Color',COL_DARK);

quiver(axA,0.29,0.73,0.075,-0.13,0,'Color',COL_GREY,'LineWidth',1.2,'MaxHeadSize',0.25);
quiver(axA,0.29,0.36,0.075,0.10,0,'Color',COL_GREY,'LineWidth',1.2,'MaxHeadSize',0.25);

rectangle(axA,'Position',[0.69 0.475 0.12 0.12],'Curvature',[1 1], ...
    'FaceColor',[1 1 1],'EdgeColor',COL_DARK,'LineWidth',1.3);
text(axA,0.75,0.535,'SV','HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontSize',10,'FontWeight','bold');
quiver(axA,0.63,0.535,0.045,0,0,'Color',COL_DARK,'LineWidth',1.2,'MaxHeadSize',0.30);

rectangle(axA,'Position',[0.86 0.63 0.11 0.18],'Curvature',0.06, ...
    'FaceColor',[0.97 0.97 0.97],'EdgeColor',COL_DARK,'LineWidth',1.1);
text(axA,0.915,0.745,'|SV|','HorizontalAlignment','center', ...
    'FontSize',9.3,'FontWeight','bold');
text(axA,0.915,0.680,{'Decision';'conflict'}, ...
    'HorizontalAlignment','center','FontSize',7.6);
quiver(axA,0.81,0.565,0.050,0.095,0,'Color',COL_DARK,'LineWidth',1.1,'MaxHeadSize',0.30);

rectangle(axA,'Position',[0.86 0.27 0.11 0.18],'Curvature',0.06, ...
    'FaceColor',[0.97 0.97 0.97],'EdgeColor',COL_GREY,'LineWidth',1.1);
text(axA,0.915,0.385,'P(accept)','HorizontalAlignment','center', ...
    'FontSize',8.6,'FontWeight','bold');
text(axA,0.915,0.315,{'Choice rule';'\tau'}, ...
    'HorizontalAlignment','center','FontSize',7.5);
quiver(axA,0.81,0.505,0.050,-0.095,0,'Color',COL_GREY,'LineWidth',1.1,'MaxHeadSize',0.30);

title(axA,'A   Computational architecture','FontSize',11.5,'FontWeight','normal');

export_one_figure(figA,out_dir,'CPT_A_ComputationalArchitecture',16.0,11.0);

figB = figure( ...
    'Name','CPT B - Probability weighting', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 12.0 10.0], ...
    'Renderer','painters');

axB = axes('Parent',figB,'Position',[0.16 0.18 0.76 0.70]);
hold(axB,'on');

plot(axB,[0 1],[0 1],'--','Color',COL_MIDGREY,'LineWidth',1.0);
plot(axB,pGrid,wGrid,'-','Color',COL_CPT,'LineWidth',2.5);

refP = [0.10 0.25 0.50 0.75 0.90];
refW = CPT_weight(refP, med_gamma);
scatter(axB,refP,refW,28,'MarkerFaceColor','w','MarkerEdgeColor',COL_CPT,'LineWidth',1.0);

xlim(axB,[0 1]);
ylim(axB,[0 1]);
pbaspect(axB,[1 1 1]);

xlabel(axB,'Objective probability','FontSize',9.2);
ylabel(axB,'Decision weight','FontSize',9.2);
title(axB,'B   Probability weighting','FontSize',11.5,'FontWeight','normal');

text(axB,0.05,0.94,sprintf('median \\gamma = %.3f',med_gamma), ...
    'Units','normalized','FontSize',8.1,'Color',COL_CPT,'VerticalAlignment','top');

set_pub_axes(axB);
loosen_axis_layout(axB, 0.10, 0.10);

export_one_figure(figB,out_dir,'CPT_B_ProbabilityWeighting',12.0,10.0);

figC = figure( ...
    'Name','CPT C - Outcome valuation', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[4 4 12.0 10.0], ...
    'Renderer','painters');

axC = axes('Parent',figC,'Position',[0.16 0.18 0.76 0.70]);
hold(axC,'on');

xline(axC,0,'-','Color',[0.82 0.82 0.82],'LineWidth',0.8);
yline(axC,0,'-','Color',[0.82 0.82 0.82],'LineWidth',0.8);

plot(axC,outcomeGrid(idxLoss),utilityGrid(idxLoss),'-','Color',COL_AMBIG,'LineWidth',2.4);
plot(axC,outcomeGrid(idxGain),utilityGrid(idxGain),'-','Color',COL_RISKY,'LineWidth',2.4);

xlabel(axC,'Objective outcome','FontSize',9.2);
ylabel(axC,'Subjective utility','FontSize',9.2);
title(axC,'C   Outcome valuation','FontSize',11.5,'FontWeight','normal');

text(axC,0.05,0.94,sprintf('\\alpha = %.3f    \\lambda = %.3f',med_alpha,med_lambda), ...
    'Units','normalized','VerticalAlignment','top','FontSize',8.1,'Color',COL_DARK);

text(axC,0.05,0.08,'Losses','Units','normalized','FontSize',8,'Color',COL_AMBIG);
text(axC,0.95,0.88,'Gains','Units','normalized','HorizontalAlignment','right','FontSize',8,'Color',COL_RISKY);

set_pub_axes(axC);
loosen_axis_layout(axC, 0.10, 0.10);

export_one_figure(figC,out_dir,'CPT_C_OutcomeValuation',12.0,10.0);

figD = figure( ...
    'Name','CPT D - Subjective-value distributions', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[5 5 16.0 10.5], ...
    'Renderer','painters');

axD = axes('Parent',figD,'Position',[0.10 0.20 0.84 0.68]);
hold(axD,'on');

risk_idx = string(T.Trial_Type) == "Risky";
amb_idx  = string(T.Trial_Type) == "Ambiguous";

svRisk = T.SV(risk_idx);
svAmb  = T.SV(amb_idx);

svRisk = svRisk(isfinite(svRisk));
svAmb  = svAmb(isfinite(svAmb));

xRisk = sign(svRisk) .* log10(1 + abs(svRisk));
xAmb  = sign(svAmb)  .* log10(1 + abs(svAmb));

xMin = min([xRisk; xAmb]);
xMax = max([xRisk; xAmb]);
xGrid = linspace(xMin, xMax, 500);

fRisk = ksdensity(xRisk, xGrid);
fAmb  = ksdensity(xAmb,  xGrid);

fRisk = fRisk ./ max(fRisk);
fAmb  = fAmb  ./ max(fAmb);

fill(axD,[xGrid fliplr(xGrid)],[fRisk zeros(size(fRisk))],COL_RISKY, ...
    'FaceAlpha',0.13,'EdgeColor','none');
fill(axD,[xGrid fliplr(xGrid)],[fAmb zeros(size(fAmb))],COL_AMBIG, ...
    'FaceAlpha',0.13,'EdgeColor','none');

plot(axD,xGrid,fRisk,'Color',COL_RISKY,'LineWidth',2.1);
plot(axD,xGrid,fAmb,'Color',COL_AMBIG,'LineWidth',2.1);

xline(axD,0,'--','Color',[0.35 0.35 0.35],'LineWidth',1.0);

rawTicks = [-100 -10 -1 0 1 10 100];
transTicks = sign(rawTicks).*log10(1+abs(rawTicks));
keep = transTicks >= xMin & transTicks <= xMax;

xticks(axD,transTicks(keep));
xticklabels(axD,string(rawTicks(keep)));
xtickangle(axD,25);

ylim(axD,[0 1.14]);
yticks(axD,[]);

xlabel(axD,'Signed subjective value (SV)','FontSize',9.5);
title(axD,'D   Empirical subjective-value distributions','FontSize',11.5,'FontWeight','normal');

text(axD,0.015,0.78,'Reject preference', ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'FontSize',7.8, ...
    'Color',[0.40 0.40 0.40], ...
    'BackgroundColor','w', ...
    'Margin',1);

text(axD,0.985,0.78,'Accept preference', ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'FontSize',7.8, ...
    'Color',[0.40 0.40 0.40], ...
    'BackgroundColor','w', ...
    'Margin',1);

hD1 = plot(axD,nan,nan,'-','Color',COL_RISKY,'LineWidth',2);
hD2 = plot(axD,nan,nan,'-','Color',COL_AMBIG,'LineWidth',2);

lgdD = legend(axD,[hD1 hD2],{'Risky','Ambiguous'}, ...
    'Box','off','FontSize',8,'Location','none');
lgdD.Units = 'normalized';
lgdD.Position = [0.02 0.79 0.20 0.12];

annotation(figD,'textbox',[0.435 0.76 0.16 0.11], ...
    'String',{'Indifference','Peak conflict'}, ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontWeight','bold', ...
    'FontSize',8.5, ...
    'BackgroundColor','w', ...
    'Margin',2);

set_pub_axes(axD);
loosen_axis_layout(axD, 0.13, 0.11);

export_one_figure(figD,out_dir,'CPT_D_SubjectiveValueDistributions',16.0,10.5);

figE = figure( ...
    'Name','CPT E - Participant-level model performance', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[6 6 16.5 10.0], ...
    'Renderer','painters');

annotation(figE,'textbox',[0.10 0.925 0.80 0.04], ...
    'String','E   Participant-level model performance', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontName','Arial', ...
    'FontSize',11.5);

% Prediction accuracy
axE1 = axes('Parent',figE,'Position',[0.10 0.18 0.36 0.63]);
hold(axE1,'on');

validAcc = isfinite(accRisky) & isfinite(accAmbig);
rAcc = accRisky(validAcc);
aAcc = accAmbig(validAcc);

for s = 1:length(rAcc)
    plot(axE1,[1 2],[rAcc(s) aAcc(s)],'-','Color',[0.86 0.86 0.86],'LineWidth',0.55);
end

rng(42,'twister');
j1 = (rand(size(rAcc))-0.5)*0.10;
j2 = (rand(size(aAcc))-0.5)*0.10;

scatter(axE1,1+j1,rAcc,20,'MarkerFaceColor',COL_RISKY_LIGHT,'MarkerEdgeColor',COL_RISKY,'LineWidth',0.45);
scatter(axE1,2+j2,aAcc,20,'MarkerFaceColor',COL_AMBIG_LIGHT,'MarkerEdgeColor',COL_AMBIG,'LineWidth',0.45);

medR = median(rAcc,'omitnan');
medA = median(aAcc,'omitnan');

plot(axE1,[0.78 1.22],[medR medR],'Color',COL_RISKY,'LineWidth',2.8);
plot(axE1,[1.78 2.22],[medA medA],'Color',COL_AMBIG,'LineWidth',2.8);

xlim(axE1,[0.55 2.45]);
ylim(axE1,[0 1.06]);
xticks(axE1,[1 2]);
xticklabels(axE1,{'Risky','Ambiguous'});
ylabel(axE1,'Prediction accuracy');
title(axE1,'Prediction accuracy','FontSize',9.5,'FontWeight','normal');

text(axE1,1,medR+0.035,sprintf('%.3f',medR), ...
    'HorizontalAlignment','center','FontSize',7.2,'FontWeight','bold','Color',COL_RISKY);
text(axE1,2,medA+0.035,sprintf('%.3f',medA), ...
    'HorizontalAlignment','center','FontSize',7.2,'FontWeight','bold','Color',COL_AMBIG);

set_pub_axes(axE1);
loosen_axis_layout(axE1, 0.09, 0.10);

% McFadden R2
axE2 = axes('Parent',figE,'Position',[0.59 0.18 0.26 0.63]);
hold(axE2,'on');

r2 = r2Risky(isfinite(r2Risky));
j3 = (rand(size(r2))-0.5)*0.16;

scatter(axE2,1+j3,r2,20,'MarkerFaceColor',COL_CPT_LIGHT,'MarkerEdgeColor',COL_CPT,'LineWidth',0.45);

medR2 = median(r2,'omitnan');
plot(axE2,[0.74 1.26],[medR2 medR2],'Color',COL_CPT,'LineWidth',2.8);

xlim(axE2,[0.55 1.45]);
ylim(axE2,[0 1.06]);
xticks(axE2,1);
xticklabels(axE2,{'McFadden R^2'});
ylabel(axE2,'Pseudo-R^2');
title(axE2,'Model fit','FontSize',9.5,'FontWeight','normal');

text(axE2,1,medR2+0.035,sprintf('%.3f',medR2), ...
    'HorizontalAlignment','center','FontSize',7.2,'FontWeight','bold','Color',COL_CPT);

set_pub_axes(axE2);
loosen_axis_layout(axE2, 0.09, 0.10);

annotation(figE,'textbox',[0.08 0.055 0.84 0.04], ...
    'String','Risky = fitted trials    |    Ambiguous = out-of-fit-set prediction', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'FontName','Arial', ...
    'FontSize',7.2, ...
    'Color',COL_GREY);

export_one_figure(figE,out_dir,'CPT_E_ParticipantModelPerformance',16.5,10.0);

figF = figure( ...
    'Name','CPT F - Choice-probability calibration', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[7 7 12.0 10.0], ...
    'Renderer','painters');

axF = axes('Parent',figF,'Position',[0.16 0.18 0.76 0.70]);
hold(axF,'on');

CalRisk = calibration_summary(T,"Risky",6);
CalAmb  = calibration_summary(T,"Ambiguous",6);

hIdent = plot(axF,[0 1],[0 1],'--','Color',[0.70 0.70 0.70],'LineWidth',1.0);

hRisk = errorbar(axF, ...
    CalRisk.Mean_Predicted, CalRisk.Mean_Observed, CalRisk.SEM_Observed, ...
    'o-','Color',COL_RISKY,'MarkerFaceColor','w','MarkerEdgeColor',COL_RISKY, ...
    'MarkerSize',5.2,'LineWidth',1.6,'CapSize',3);

hAmb = errorbar(axF, ...
    CalAmb.Mean_Predicted, CalAmb.Mean_Observed, CalAmb.SEM_Observed, ...
    'o-','Color',COL_AMBIG,'MarkerFaceColor','w','MarkerEdgeColor',COL_AMBIG, ...
    'MarkerSize',5.2,'LineWidth',1.6,'CapSize',3);

xlim(axF,[0 1]);
ylim(axF,[0 1]);
pbaspect(axF,[1 1 1]);

xlabel(axF,'Model-predicted P(accept)','FontSize',9.2);
ylabel(axF,'Observed P(accept)','FontSize',9.2);
title(axF,'F   Choice-probability calibration','FontSize',11.5,'FontWeight','normal');

legend(axF,[hIdent hRisk hAmb],{'Identity','Risky','Ambiguous'}, ...
    'Location','northwest','Box','off','FontSize',8);

set_pub_axes(axF);
loosen_axis_layout(axF, 0.10, 0.10);

export_one_figure(figF,out_dir,'CPT_F_ChoiceProbabilityCalibration',12.0,10.0);

ProbabilitySource = table(pGrid, wGrid, ...
    'VariableNames', {'ObjectiveProbability','DecisionWeight'});
writetable(ProbabilitySource, fullfile(out_dir,'CPT_B_ProbabilityWeighting_Source.csv'));

ValueSource = table(outcomeGrid, utilityGrid, ...
    'VariableNames', {'ObjectiveOutcome','SubjectiveUtility'});
writetable(ValueSource, fullfile(out_dir,'CPT_C_OutcomeValuation_Source.csv'));

fprintf('Figures exported\n');
fprintf('Output folder:\n%s\n', out_dir);
fprintf('Formats: FIG, PDF, EPS, PNG, TIF');
fprintf(' and SVG if exportgraphics is supported.\n');

function w = CPT_weight(p, gamma)
    p = max(min(p, 1-1e-12), 1e-12);
    w = p.^gamma ./ (p.^gamma + (1-p).^gamma).^(1/gamma);
end

function set_pub_axes(ax)
    set(ax, ...
        'FontName','Arial', ...
        'FontSize',8.4, ...
        'LineWidth',0.85, ...
        'TickDir','out', ...
        'TickLength',[0.018 0.018], ...
        'Box','off', ...
        'Layer','top');
end

function loosen_axis_layout(ax, xShiftDown, yShiftLeft)
    % Adjust axis spacing.
    pos = ax.Position;
    pos(1) = pos(1) + 0.01;
    pos(2) = pos(2) + 0.02;
    pos(3) = pos(3) - 0.02;
    pos(4) = pos(4) - 0.02;
    ax.Position = pos;

    drawnow;

    % Move x-label.
    ax.XLabel.Units = 'normalized';
    xl = ax.XLabel.Position;
    xl(2) = xl(2) - xShiftDown;
    ax.XLabel.Position = xl;

    % Move y-label.
    ax.YLabel.Units = 'normalized';
    yl = ax.YLabel.Position;
    yl(1) = yl(1) - yShiftLeft;
    ax.YLabel.Position = yl;
end

function export_one_figure(fig, out_dir, baseName, w_cm, h_cm)

    drawnow;
    pause(0.20);

    % Save MATLAB figure.
    savefig(fig, fullfile(out_dir,[baseName '.fig']));

    % Set paper size.
    set(fig,'PaperUnits','centimeters');
    set(fig,'PaperSize',[w_cm h_cm]);
    set(fig,'PaperPosition',[0 0 w_cm h_cm]);
    set(fig,'PaperPositionMode','manual');
    set(fig,'InvertHardcopy','off');

    pdf_file = fullfile(out_dir,[baseName '.pdf']);
    eps_file = fullfile(out_dir,[baseName '.eps']);
    png_file = fullfile(out_dir,[baseName '.png']);
    tif_file = fullfile(out_dir,[baseName '.tif']);

    print(fig,pdf_file,'-dpdf','-painters');
    print(fig,eps_file,'-depsc','-painters');
    print(fig,png_file,'-dpng','-r600');
    print(fig,tif_file,'-dtiff','-r600');

    % Export SVG when available.
    try
        exportgraphics(fig, fullfile(out_dir,[baseName '.svg']), ...
            'ContentType','vector');
    catch
        % SVG export is optional.
    end
end

function Summary = calibration_summary(T, conditionName, nBins)

    conditionMask = ...
        string(T.Trial_Type) == conditionName & ...
        isfinite(T.P_Accept) & ...
        isfinite(T.Choice);

    D = T(conditionMask,:);
    p = D.P_Accept;

    edges = quantile(p, linspace(0,1,nBins+1));
    edges = unique(edges, 'stable');

    edges(1)   = min(p) - eps;
    edges(end) = max(p) + eps;

    if length(edges) < 3
        error('Too few unique predicted probabilities for calibration.');
    end

    bin = discretize(p, edges);
    IDs = unique(string(D.ID));

    OUT_ID   = strings(0,1);
    OUT_BIN  = zeros(0,1);
    OUT_PRED = zeros(0,1);
    OUT_OBS  = zeros(0,1);
    OUT_N    = zeros(0,1);

    for s = 1:length(IDs)
        sid = IDs(s);
        idxS = string(D.ID) == sid;
        binsPresent = unique(bin(idxS & ~isnan(bin)));

        for b = binsPresent'
            idx = idxS & bin == b;

            if sum(idx) < 1
                continue;
            end

            OUT_ID(end+1,1)   = sid;
            OUT_BIN(end+1,1)  = b;
            OUT_PRED(end+1,1) = mean(D.P_Accept(idx));
            OUT_OBS(end+1,1)  = mean(D.Choice(idx));
            OUT_N(end+1,1)    = sum(idx);
        end
    end

    ParticipantBins = table( ...
        OUT_ID, OUT_BIN, OUT_PRED, OUT_OBS, OUT_N, ...
        'VariableNames', {'ID','Bin','Predicted','Observed','NTrials'});

    validBins = unique(ParticipantBins.Bin);

    BIN    = zeros(length(validBins),1);
    MEAN_P = zeros(length(validBins),1);
    MEAN_O = zeros(length(validBins),1);
    SEM_O  = zeros(length(validBins),1);
    NSUB   = zeros(length(validBins),1);

    for k = 1:length(validBins)
        b = validBins(k);
        idx = ParticipantBins.Bin == b;

        BIN(k)    = b;
        MEAN_P(k) = mean(ParticipantBins.Predicted(idx));
        MEAN_O(k) = mean(ParticipantBins.Observed(idx));
        NSUB(k)   = sum(idx);

        SEM_O(k)  = std(ParticipantBins.Observed(idx)) / sqrt(NSUB(k));
    end

    Summary = table( ...
        BIN, MEAN_P, MEAN_O, SEM_O, NSUB, ...
        'VariableNames', {'Bin','Mean_Predicted','Mean_Observed','SEM_Observed','N_Subjects'});
end
