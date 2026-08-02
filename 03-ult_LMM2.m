
% LMM2: THE NEUROCOMPUTATIONAL BRIDGE (PART 3)

clear; close all;
%%
disp('=================================================================');
disp('STAGE 1: DATA PREP & GLOBAL Z-SCORING');
disp('=================================================================');
% 1. Load Behavioral Data
behav_data = readtable('analysis_ready_data_extreme_groups.csv'); % Assuming this has the CPT parameters
behav_data = behav_data(~isnan(behav_data.RT) & ~isnan(behav_data.Abs_SV), :);


% --- CRITICAL FIX: FORCING THE CORRECT 19 vs 19 LABELS ---

high_ids = [6, 17, 21, 22, 25, 29, 31, 32, 33, 36, 37, 39, 41, 46, 49, 54, 55, 58, 61];
low_ids  = [1, 9, 11, 12, 13, 15, 27, 28, 38, 40, 45, 47, 48, 50, 52, 53, 60, 66, 67];

% Convert CC_Group to string temporarily to easily overwrite them
behav_data.CC_Group = string(behav_data.CC_Group);

% Reset EVERYONE to 'Med CC' first (to clear old labels)
behav_data.CC_Group(:) = "Med CC";

% Apply the true labels based on your verified correct list
behav_data.CC_Group(ismember(double(string(behav_data.ID)), low_ids)) = "Low CC";
behav_data.CC_Group(ismember(double(string(behav_data.ID)), high_ids)) = "High CC";

% Convert back to categorical
behav_data.CC_Group = categorical(behav_data.CC_Group);
% =========================================================================

% 2. Load EEG ROI Data
opts = detectImportOptions('ROI_4Cond_pt_beta_all2.xlsx');
opts.VariableNamingRule = 'preserve'; 
eeg_data = readtable('ROI_4Cond_pt_beta_all2.xlsx', opts);
% 3. Calculate Trait FPN Neural Braking (HA minus RP)
FPN_RP = mean(eeg_data{:, {'P3_RP', 'P4_RP', 'Pz_RP', 'CPz_RP', 'POz_RP', 'F3_RP', 'F4_RP', 'Fz_RP', 'FCz_RP'}}, 2);
FPN_HA = mean(eeg_data{:, {'P3_HA', 'P4_HA', 'Pz_HA', 'CPz_HA', 'POz_HA', 'F3_HA', 'F4_HA', 'Fz_HA', 'FCz_HA'}}, 2);
eeg_data.Neural_Braking = FPN_HA - FPN_RP;
% 4. Inject Neural Braking into Behavioral Data
behav_data.Neural_Braking = NaN(height(behav_data), 1);
for s = 1:height(eeg_data)
    idx = (behav_data.ID == eeg_data.Subject_ID(s));
    behav_data.Neural_Braking(idx) = eeg_data.Neural_Braking(s);
end
% 5. Isolate Ambiguous Trials ONLY
ambig_data = behav_data(strcmp(behav_data.Trial_Type, 'Ambiguous'), :);
ambig_data = ambig_data(~isnan(ambig_data.Neural_Braking), :);
% --- CRITICAL FIX 1: GLOBAL Z-SCORING BEFORE FILTERING EXTREMES ---
ambig_data.Log_RT = log(ambig_data.RT);
ambig_data.Z_Abs_SV = (ambig_data.Abs_SV - mean(ambig_data.Abs_SV)) / std(ambig_data.Abs_SV);
ambig_data.Z_Neural_Braking = (ambig_data.Neural_Braking - mean(ambig_data.Neural_Braking)) / std(ambig_data.Neural_Braking);

% 6. Drop the Med CC Group for Extreme Contrast
extreme_data = ambig_data(ambig_data.CC_Group == 'High CC' | ambig_data.CC_Group == 'Low CC', :);
extreme_data.CC_Group = removecats(extreme_data.CC_Group, 'Med CC'); 
extreme_data.CC_Group = reordercats(extreme_data.CC_Group, {'Low CC', 'High CC'});
disp('Data prep complete. Global Z-scoring applied successfully.');
disp('=================================================================');
disp('STAGE 2: THE 3-WAY NEUROCOMPUTATIONAL LMM (HIGH vs. LOW)');
disp('=================================================================');
disp('Formula: Log_RT ~ Z_Abs_SV * Z_Neural_Braking * CC_Group + (Z_Abs_SV | ID)');
% Fit Full 3-Way Model
lmm_extreme = fitlme(extreme_data, ...
    'Log_RT ~ Z_Abs_SV * Z_Neural_Braking * CC_Group + (Z_Abs_SV | ID)', ...
    'FitMethod', 'REML');
% Print the ANOVA table
disp('--- FULL 3-WAY ANOVA RESULTS ---');
disp(anova(lmm_extreme));
% --- CRITICAL FIX 2: POST-HOC SIMPLE SLOPES ANALYSIS
disp('=================================================================');
disp('STAGE 2.5: SIMPLE SLOPES DECOMPOSITION (THE DECOUPLING TEST)');
disp('=================================================================');
% Test High CC Group
high_data = extreme_data(extreme_data.CC_Group == 'High CC', :);
lmm_high = fitlme(high_data, 'Log_RT ~ Z_Abs_SV * Z_Neural_Braking + (Z_Abs_SV | ID)', 'FitMethod', 'REML');
disp('--- HIGH CC GROUP: NEUROCOMPUTATIONAL BRIDGE ---');
disp(anova(lmm_high));
% Test Low CC Group
low_data = extreme_data(extreme_data.CC_Group == 'Low CC', :);
lmm_low = fitlme(low_data, 'Log_RT ~ Z_Abs_SV * Z_Neural_Braking + (Z_Abs_SV | ID)', 'FitMethod', 'REML');
disp('--- LOW CC GROUP: DECOUPLED SYSTEM ---');
disp(anova(lmm_low));
% -----------------------------------------------------------------------
disp('=================================================================');
disp('STAGE 3: CALCULATING MODEL ACCURACY (R-SQUARED)');
disp('=================================================================');
y_actual = extreme_data.Log_RT;
SST = sum((y_actual - mean(y_actual)).^2);
% Conditional R-squared (Fixed + Random Effects)
y_pred_cond = predict(lmm_extreme); 
SSR_cond = sum((y_actual - y_pred_cond).^2);
R2_Conditional = 1 - (SSR_cond / SST);
% Marginal R-squared (Fixed Effects ONLY)
y_pred_marg = predict(lmm_extreme, 'Conditional', false);
SSR_marg = sum((y_actual - y_pred_marg).^2);
R2_Marginal = 1 - (SSR_marg / SST);
fprintf('Marginal R-squared (Theory Only): %.2f%%\n', R2_Marginal * 100);
fprintf('Conditional R-squared (Theory + Subjects): %.2f%%\n', R2_Conditional * 100);
disp('=================================================================');
disp('STAGE 4: PUBLICATION-READY 2D INTERACTION PLOT');
disp('=================================================================');
figure('Color', 'w', 'Position', [100, 100, 1050, 500], 'Name', 'Neurocomputational Bridge');
groups = {'Low CC', 'High CC'};
beta_levels = [1, -1]; 
beta_labels = {'High Beta Capacity (+1 SD)', 'Low Beta Capacity (-1 SD)'};
colors = {[0.85, 0.25, 0.15], [0.15, 0.45, 0.75]}; 
x_sv_z = linspace(-2, 2, 100)'; 
for g = 1:2
    subplot(1, 2, g);
    hold on; grid on;
    title(sprintf('%s Group', groups{g}), 'FontSize', 15);
    
    for b = 1:2
        T_pred = table(...
            repmat(categorical(groups(g)), 100, 1), x_sv_z, repmat(beta_levels(b), 100, 1), ...
            repmat(extreme_data.ID(1), 100, 1), ... 
            'VariableNames', {'CC_Group', 'Z_Abs_SV', 'Z_Neural_Braking', 'ID'});
        
        y_pred = predict(lmm_extreme, T_pred, 'Conditional', false);
        plot(x_sv_z, y_pred, 'Color', colors{b}, 'LineWidth', 3.5, 'DisplayName', beta_labels{b});
    end
    
    xlabel('Math Conflict (Z-Scored SV)', 'FontWeight', 'bold');
    ylabel('Log Reaction Time (s)', 'FontWeight', 'bold');
    
   
    ylim([-2, 2]); 
    
    set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'off', 'TickDir', 'out');
    if g == 2
        legend('Location', 'southeast', 'Box', 'off', 'FontSize', 11);
    end
end
sgtitle('Neurocomputational Bridge: FPN Beta Capacity Moderates Decision Threshold', 'FontSize', 16, 'FontWeight', 'bold');
disp('=================================================================');
disp('STAGE 5: PUBLICATION-READY 3D SURFACE PLOT');
disp('=================================================================');
figure('Color', 'w', 'Position', [50, 50, 1200, 550], 'Name', '3D Neurocomputational Bridge');
x_sv_vec = linspace(-2, 2, 30); 
y_beta_vec = linspace(-2, 2, 30); 
[X_SV, Y_Beta] = meshgrid(x_sv_vec, y_beta_vec);
pred_SV = X_SV(:); pred_Beta = Y_Beta(:);
cmap_low = flipud(autumn(256)); 
cmap_high = winter(256);        
for g = 1:2
    subplot(1, 2, g);
    
    T_pred = table(...
        repmat(categorical(groups(g)), length(pred_SV), 1), pred_SV, pred_Beta, ...
        repmat(extreme_data.ID(1), length(pred_SV), 1), ... 
        'VariableNames', {'CC_Group', 'Z_Abs_SV', 'Z_Neural_Braking', 'ID'});
    
    Z_RT_pred = predict(lmm_extreme, T_pred, 'Conditional', false);
    Z_RT_Matrix = reshape(Z_RT_pred, size(X_SV));
    
    surf(X_SV, Y_Beta, Z_RT_Matrix, 'EdgeAlpha', 0.4, 'FaceAlpha', 0.85);
    hold on;
    
    if g == 1; colormap(gca, cmap_low); else; colormap(gca, cmap_high); end
    
    cb = colorbar; cb.Label.String = 'Log Reaction Time'; cb.Label.FontWeight = 'bold';
    
    zlim([-2, 2]); caxis([-2, 2]); % Adjusted to match 2D plot limits
    
    xlabel('Math Conflict (Z-SV)', 'FontWeight', 'bold', 'Rotation', 15);
    ylabel('Beta Capacity (Z)', 'FontWeight', 'bold', 'Rotation', -25);
    zlabel('Log Reaction Time (s)', 'FontWeight', 'bold');
    title(sprintf('%s Group', groups{g}), 'FontSize', 16);
    
    view(-45, 25); grid on;
    set(gca, 'FontSize', 11, 'LineWidth', 1);
    camlight('headlight'); lighting gouraud;
end
sgtitle('3D Neurocomputational Space: FPN Beta Moderates the Decision Threshold', 'FontSize', 18, 'FontWeight', 'bold');
rotate3d on;
%%
% =========================================================================
% STAGE 4:  2D INTERACTION (WITH 95% CIs)
% =========================================================================
disp('=================================================================');
disp('GENERATING Q1 PUBLICATION FIGURE (2D + CONFIDENCE INTERVALS)');
disp('=================================================================');
if ~exist('extreme_data', 'var') || ~exist('lmm_extreme', 'var')
    error('Please run the main LMM analysis first.');
end

% 1. Some SETTINGS

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultTextFontName', 'Arial');
% Classic colorblind-friendly scientific palette
color_high_beta = [0.85, 0.32, 0.09];  % Deep Orange-Red
color_low_beta  = [0.00, 0.44, 0.74];  % Deep Blue
groups = {'Low CC', 'High CC'};
beta_levels = [1, -1]; 
beta_labels = {'High β Capacity (+1 SD)', 'Low β Capacity (−1 SD)'};
% X-Axis: 100 points for a smooth line
x_sv_z = linspace(-2, 2, 100)';
% True Dynamic Y-Axis Limits from ACTUAL Data
min_log_rt = min(extreme_data.Log_RT) - 0.2;
max_log_rt = max(extreme_data.Log_RT) + 0.2;

% 2. GENERATE THE FIGURE
% ---------------------------------------------------------
fig1 = figure('Color', 'w', 'Position', [100, 100, 900, 450], 'Name', 'Neurocomputational Bridge');
for g = 1:2
    subplot(1, 2, g);
    hold on; box off;
    
    % Arrays to store handles for the legend
    line_handles = zeros(1, 2);
    
    for b = 1:2
        % Create prediction table
        T_pred = table(...
            repmat(categorical(groups(g)), 100, 1), x_sv_z, repmat(beta_levels(b), 100, 1), ...
            repmat(extreme_data.ID(1), 100, 1), ...
            'VariableNames', {'CC_Group', 'Z_Abs_SV', 'Z_Neural_Braking', 'ID'});
        
        % --- Extract predictions and 95% confidence intervals ---
        [y_pred, y_ci] = predict(lmm_extreme, T_pred, 'Conditional', false);
        
        % Assign colors
        if b == 1
            line_color = color_high_beta;
        else
            line_color = color_low_beta;
        end
        
        % 1. Plot the Shaded 95% Confidence Interval Ribbon
        x_fill = [x_sv_z; flipud(x_sv_z)];
        y_fill = [y_ci(:, 1); flipud(y_ci(:, 2))]; % Lower bound then upper bound reversed
        
        fill(x_fill, y_fill, line_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off'); % Don't show the ribbon box in the legend
         
        % 2. Plot the Main Prediction Line
        line_handles(b) = plot(x_sv_z, y_pred, '-', 'Color', line_color, 'LineWidth', 3.5);
    end
    
   
    % Axes Labels
    xlabel('Mathematical Conflict (Standardized Z-|SV|)', 'FontWeight', 'bold');
    if g == 1
        ylabel('Log Reaction Time', 'FontWeight', 'bold');
    end
    
    % Title
    title(sprintf('%s Phenotype', groups{g}), 'FontSize', 15, 'FontWeight', 'bold');
    
    % Semantic Text: Explicitly defining what -2 and +2 mean for the Reviewer
    yl = [min_log_rt, max_log_rt];
    text(-1.9, yl(2)*0.98, 'Peak Conflict', 'FontSize', 10, 'Color', [0.4 0.4 0.4]);
    text(1.9, yl(2)*0.98, 'No Conflict', 'FontSize', 10, 'Color', [0.4 0.4 0.4], 'HorizontalAlignment', 'right');
    
    % Apply Limits
    ylim(yl);
    xlim([-2, 2]);
    
    % Clean Grid and Ticks
    grid on; 
    set(gca, 'GridAlpha', 0.15, 'TickDir', 'out', 'LineWidth', 1.2);
    
    % Legend (Only on the right panel)
    if g == 2
        lgd = legend(line_handles, beta_labels, 'Location', 'northeast', 'Box', 'off', 'FontSize', 11);
    end
end
% Add a clean title
sgtitle('Neurocomputational Decoupling: FPN Beta Capacity Moderates Decision Threshold', ...
        'FontSize', 16, 'FontWeight', 'bold');
% Save High-resolution outputs
exportgraphics(fig1, 'Fig1_Interaction_95CI.png', 'Resolution', 600);
exportgraphics(fig1, 'Fig1_Interaction_95CI.pdf', 'ContentType', 'vector'); % PDFs are required by many Q1 journals
fprintf('   ✅ Figure saved successfully (PNG and PDF formats).\n');
