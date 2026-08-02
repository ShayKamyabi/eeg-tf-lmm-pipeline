
%  LMM ANALYSIS: THE TEMPORAL BRAKING MECHANISM 

clear; close all; clc;

%% 1. Load Data
disp('=================================================================');
disp('STAGE 1: DATA PREP & FORCING CORRECT 19 vs 19 LABELS');
disp('=================================================================');
disp('Loading master_data_fixed_v4.csv...');
data = readtable('master_data_fixed_v4.csv');

% Clean data
data = data(~isnan(data.RT), :);

% Log-Transform RT for Statistical Consistency
data.Log_RT = log(data.RT);
data.ID = categorical(data.ID);


% --- CRITICAL FIX: FORCING THE CORRECT 19 vs 19 LABELS ---

% Define the exact IDs based on your updated Stroop boundaries
low_ids  = [1, 9, 11, 12, 13, 15, 27, 28, 38, 40, 45, 47, 48, 50, 52, 53, 60, 66, 67];
high_ids = [6, 17, 21, 22, 25, 29, 31, 32, 33, 36, 37, 39, 41, 46, 49, 54, 55, 58, 61];

% Convert CC_Group to string temporarily to easily overwrite them
data.CC_Group = string(data.CC_Group);

% Reset EVERYONE to 'Med CC' first (to clear old wrong labels)
data.CC_Group(:) = "Med CC";

% Apply the true labels
data.CC_Group(ismember(double(string(data.ID)), true_low_ids)) = "Low CC";
data.CC_Group(ismember(double(string(data.ID)), true_high_ids)) = "High CC";

% Convert back to categorical 
data.CC_Group = categorical(data.CC_Group);
% =========================================================================

%% 2. COMPUTATIONAL SV EXTRACTION (DONE GLOBALLY FIRST)
disp('=================================================================');
disp('STAGE 2: EXTRACTING CPT PARAMETERS & SV CALCULATIONS');
disp('=================================================================');
subject_ids = unique(data.ID);
data.SV = NaN(height(data), 1);
options = optimoptions('fmincon', 'Display', 'off');
lb = [0.01, 0.01, 0.01, 0.01]; ub = [2.00, 10.00, 2.00, 10.00];
init_guess = [0.88, 2.25, 0.65, 1.0];

for s = 1:length(subject_ids)
    current_id = subject_ids(s);
    idx_sub = data.ID == current_id;
    
    % Fit parameters on Risky trials only to get a "pure" baseline
    idx_fit = idx_sub & strcmp(data.Trial_Type, 'Risky');
    sub_fit = data(idx_fit, :);
    
    if height(sub_fit) > 10
        obj_func = @(params) calc_NLL(params, sub_fit.Loss_Mag, sub_fit.Loss_Prob, ...
                                      sub_fit.Gain_Prob, sub_fit.Gain_Mag, sub_fit.Choice);
        best_params = fmincon(obj_func, init_guess, [], [], [], [], lb, ub, [], options);
        
        alpha = best_params(1); lambda = best_params(2); gamma = best_params(3);
        
        % Calculate SV for ALL trials
        sub_all = data(idx_sub, :);
        w_g = (sub_all.Gain_Prob.^gamma) ./ ((sub_all.Gain_Prob.^gamma + (1 - sub_all.Gain_Prob).^gamma).^(1/gamma));
        w_l = (sub_all.Loss_Prob.^gamma) ./ ((sub_all.Loss_Prob.^gamma + (1 - sub_all.Loss_Prob).^gamma).^(1/gamma));
        data.SV(idx_sub) = w_g .* (sub_all.Gain_Mag.^alpha) - lambda .* w_l .* (sub_all.Loss_Mag.^alpha);
    end
end

% The Independent Variable: Absolute SV
data.Abs_SV = abs(data.SV);
data = data(~isnan(data.Abs_SV), :);

% --- CRITICAL FIX 1: Globally Z-score |SV| before removing Med CC ---
data.Z_Abs_SV = (data.Abs_SV - mean(data.Abs_SV)) / std(data.Abs_SV);


% filter out the Med CC rows for the Extreme Groups analysis
data = data(data.CC_Group ~= 'Med CC', :);
data.CC_Group = removecats(data.CC_Group, 'Med CC'); 
data.CC_Group = reordercats(data.CC_Group, {'Low CC', 'High CC'});

%  SUBJECT COUNT VALIDATION
unique_high_subjects = length(unique(data.ID(data.CC_Group == 'High CC')));
unique_low_subjects  = length(unique(data.ID(data.CC_Group == 'Low CC')));

fprintf('\n SAMPLE SIZE VALIDATION:\n');
fprintf('   - High CC Subjects Found: %d\n', unique_high_subjects);
fprintf('   - Low CC Subjects Found:  %d\n', unique_low_subjects);
if unique_high_subjects ~= 19 || unique_low_subjects ~= 19
    fprintf(' WARNING: Sample size is not 19 vs 19.\n');
else
    disp('   - Perfect! Dataset correctly merged for 38 subjects.');
end

% Save the calculated SV data 
writetable(data, 'analysis_ready_data_extreme_groups.csv');
disp('analysis-ready dataset saved successfully!');

%% 3. FINAL LINEAR MIXED-EFFECTS MODELS (WITH RANDOM SLOPES)
disp('=================================================================');
disp('STAGE 3: BEHAVIORAL LMM MODELS (HIGH vs. LOW CC)');
disp('=================================================================');

% SPLIT 1: RISKY TRIALS
disp('--- LMM: RISKY TRIALS (Known Probabilities) ---');
risky_data = data(strcmp(data.Trial_Type, 'Risky'), :);
lmm_risky = fitlme(risky_data, 'Log_RT ~ Z_Abs_SV * CC_Group + (Z_Abs_SV | ID)', 'FitMethod', 'REML');
disp(anova(lmm_risky))

% SPLIT 2: AMBIGUOUS TRIALS 
disp('--- LMM: AMBIGUOUS TRIALS (Missing Probabilities) ---');
ambig_data = data(strcmp(data.Trial_Type, 'Ambiguous'), :);
lmm_ambig = fitlme(ambig_data, 'Log_RT ~ Z_Abs_SV * CC_Group + (Z_Abs_SV | ID)', 'FitMethod', 'REML');
disp(anova(lmm_ambig))

%% 4. PUBLICATION-READY PLOT
disp('=================================================================');
disp('GENERATING FIGURE (2D + CONFIDENCE INTERVALS)');
disp('=================================================================');

% Strict Publication Settings
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 12);

fig2 = figure('Color', 'w', 'Position', [100, 100, 1100, 500], 'Name', 'Behavioral Divergence');
subplot(1,2,1); plot_lmm(risky_data, lmm_risky, 'Risky Condition (Known Probabilities)');
subplot(1,2,2); plot_lmm(ambig_data, lmm_ambig, 'Ambiguous Condition (Peak Conflict)');
sgtitle('Cognitive Control Moderates Behavioral Braking Under Uncertainty', 'FontSize', 16, 'FontWeight', 'bold');

% Save Figure
exportgraphics(fig2, 'Fig2_Behavioral_Interaction.png', 'Resolution', 600);
exportgraphics(fig2, 'Fig2_Behavioral_Interaction.pdf', 'ContentType', 'vector');
fprintf(' Figure saved successfully (PNG and PDF formats).\n');

%% =========================================================================
% HELPER FUNCTIONS

function nll = calc_NLL(p, lm, lp, gp, gm, ch)
    w_g = (gp.^p(3))./((gp.^p(3)+(1-gp).^p(3)).^(1/p(3)));
    w_l = (lp.^p(3))./((lp.^p(3)+(1-lp).^p(3)).^(1/p(3)));
    sv = w_g.*(gm.^p(1)) - p(2).*w_l.*(lm.^p(1));
    pr = 1./(1+exp(-sv./p(4)));
    pr = max(min(pr, 1-1e-10), 1e-10);
    nll = -sum(ch.*log(pr) + (1-ch).*log(1-pr));
end

function plot_lmm(sub_data, model, titlestr)
    hold on; box off;
    
    %  Color Palette
    colors = {[0.00, 0.44, 0.74], [0.85, 0.32, 0.09]}; % Blue (Low CC), Orange (High CC)
    groups = {'Low CC', 'High CC'};
    
    x_ax = linspace(-2, 2, 100)';
    
    % Dynamic Y-Limits
    min_log_rt = min(sub_data.Log_RT) - 0.2;
    max_log_rt = max(sub_data.Log_RT) + 0.2;
    
    line_handles = zeros(1, 2);
    
    for i = 1:2
        % Generate fixed-effect predictions and Confidence Intervals
        T = table(repmat(categorical(groups(i)), 100, 1), x_ax, repmat(sub_data.ID(1), 100, 1), ...
            'VariableNames', {'CC_Group', 'Z_Abs_SV', 'ID'});
        
        [ypred, yci] = predict(model, T, 'Conditional', false);
        
        % Plot 95% CI Ribbons
        x_fill = [x_ax; flipud(x_ax)];
        y_fill = [yci(:,1); flipud(yci(:,2))];
        fill(x_fill, y_fill, colors{i}, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % Plot Main Line
        line_handles(i) = plot(x_ax, ypred, 'Color', colors{i}, 'LineWidth', 3.5, 'DisplayName', groups{i});
    end
    
    xlabel('Mathematical Conflict (Standardized Z-|SV|)', 'FontWeight', 'bold'); 
    ylabel('Log Reaction Time', 'FontWeight', 'bold'); 
    title(titlestr, 'FontSize', 13, 'FontWeight', 'bold'); 
    
    % Semantic Annotations 
    yl = [min_log_rt, max_log_rt];
    ylim(yl); xlim([-2, 2]);
    text(-1.9, yl(2)*0.98, 'Peak Conflict', 'FontSize', 10, 'Color', [0.4 0.4 0.4]);
    text(1.9, yl(2)*0.98, 'No Conflict', 'FontSize', 10, 'Color', [0.4 0.4 0.4], 'HorizontalAlignment', 'right');
    
    set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'GridAlpha', 0.15);
    grid on;
    
    legend(line_handles, groups, 'Location', 'northeast', 'Box', 'off', 'FontSize', 11);
end
