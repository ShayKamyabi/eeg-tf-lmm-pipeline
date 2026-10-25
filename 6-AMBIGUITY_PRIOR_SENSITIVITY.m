clear;
clc;
close all;

behavior_file = 'master_data_fixed_v4.csv';

% ROI9 summary used for LMM2.
roi_file = 'F:\my data\Figures\TF_MASTER_SIMPLE_FINAL\ROI9_ROI6_summary_FINAL.csv';

out_dir = 'F:\my data\Figures\Ambiguity_Prior_Sensitivity';

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

if ~isfile(behavior_file)
    [f,p] = uigetfile('*.csv','Select master_data_fixed_v4.csv');
    if isequal(f,0), error('Behavior file not selected.'); end
    behavior_file = fullfile(p,f);
end

if ~isfile(roi_file)
    [f,p] = uigetfile('*.csv','Select ROI_4Cond_pt_beta_all2.csv');
    if isequal(f,0), error('ROI file not selected.'); end
    roi_file = fullfile(p,f);
end

LOW_IDS = [ ...
     1  2 11 12 13 15 27 28 38 40 ...
    45 47 48 50 52 53 60 66 67];

HIGH_IDS = [ ...
     6 17 21 22 25 29 31 32 33 36 ...
    37 39 41 46 49 54 55 58 61];

EXTREME_IDS = sort([LOW_IDS HIGH_IDS]);

B = readtable(behavior_file);

requiredVars = {'ID','Trial_Type','Choice','RT', ...
    'Gain_Mag','Gain_Prob','Loss_Mag','Loss_Prob'};

for k = 1:numel(requiredVars)
    if ~ismember(requiredVars{k},B.Properties.VariableNames)
        error('Missing variable: %s',requiredVars{k});
    end
end

B.ID = double(B.ID);
B.Trial_Type = strtrim(string(B.Trial_Type));

% Remove ID 9.
B(B.ID==9,:) = [];

if numel(unique(B.ID)) ~= 57
    error('Expected N57 after removing ID9; found N=%d.',numel(unique(B.ID)));
end

% Assign final CC groups.
grp = repmat("Med CC",height(B),1);
grp(ismember(B.ID,LOW_IDS)) = "Low CC";
grp(ismember(B.ID,HIGH_IDS)) = "High CC";
B.CC_Group_FINAL = grp;

% Apply RT corrections.
B.RT_Original = B.RT;

rt_fix_mask = isfinite(B.RT) & B.RT > 0 & B.RT < 1;
fprintf('Sub-1-s RT rows corrected: %d\n',sum(rt_fix_mask));
if sum(rt_fix_mask) ~= 8
    warning('Expected 8 sub-1-s RT rows; found %d.',sum(rt_fix_mask));
end
B.RT(rt_fix_mask) = 10 .* B.RT(rt_fix_mask);

bad_rt_mask = B.ID==67 & strcmpi(B.Trial_Type,'Ambiguous') & ...
    isfinite(B.RT_Original) & ...
    abs(B.RT_Original - 1.00135137019606) < 1e-10;

if sum(bad_rt_mask) ~= 1
    error('Expected exactly one confirmed bad ID67 ambiguous RT row; found %d.', ...
        sum(bad_rt_mask));
end
B(bad_rt_mask,:) = [];

B.Log_RT = log(B.RT);

fprintf('Final behavior pool: N=%d, rows=%d\n',numel(unique(B.ID)),height(B));

isAmb = strcmpi(B.Trial_Type,'Ambiguous');
isRisk = strcmpi(B.Trial_Type,'Risky');

fprintf('Ambiguity input check\n');

fprintf('Unique ambiguous Gain_Prob:\n');
disp(unique(B.Gain_Prob(isAmb))');

fprintf('Unique ambiguous Loss_Prob:\n');
disp(unique(B.Loss_Prob(isAmb))');

fprintf('Ambiguous Gain_Mag == 12.5 rows: %d\n', ...
    sum(isAmb & abs(B.Gain_Mag-12.5)<1e-12));
fprintf('Ambiguous Loss_Mag == 12.5 rows: %d\n', ...
    sum(isAmb & abs(B.Loss_Mag-12.5)<1e-12));

if any(abs(B.Gain_Mag(isRisk)-12.5)<1e-12) || ...
   any(abs(B.Loss_Mag(isRisk)-12.5)<1e-12)
    warning('12.5 also appears in risky trials; verify before treating it as an imputed hidden magnitude.');
end

subject_ids = sort(unique(B.ID));

options = optimoptions('fmincon','Display','off');
lb = [0.01, 0.01, 0.01, 0.01];
ub = [2.00,10.00,2.00,10.00];
init_guess = [0.88,2.25,0.65,1.0];

Params = table(subject_ids, ...
    nan(numel(subject_ids),1), ...
    nan(numel(subject_ids),1), ...
    nan(numel(subject_ids),1), ...
    nan(numel(subject_ids),1), ...
    'VariableNames',{'ID','alpha','lambda','gamma','tau'});

for s = 1:numel(subject_ids)

    sid = subject_ids(s);

    sub_fit = B(B.ID==sid & isRisk,:);

    if height(sub_fit) < 11
        error('ID %d has fewer than 11 risky trials.',sid);
    end

    obj_func = @(p) calc_NLL( ...
        p, ...
        sub_fit.Loss_Mag, ...
        sub_fit.Loss_Prob, ...
        sub_fit.Gain_Prob, ...
        sub_fit.Gain_Mag, ...
        sub_fit.Choice);

    best_params = fmincon( ...
        obj_func,init_guess,[],[],[],[],lb,ub,[],options);

    Params.alpha(s)  = best_params(1);
    Params.lambda(s) = best_params(2);
    Params.gamma(s)  = best_params(3);
    Params.tau(s)    = best_params(4);

end

fprintf('CPT fitted for %d participants.\n',height(Params));

opts = detectImportOptions(roi_file);
opts.VariableNamingRule = 'preserve';
E = readtable(roi_file,opts);

if ~ismember('Subject_ID',E.Properties.VariableNames)
    error('ROI file is missing Subject_ID.');
end

E.Subject_ID = double(E.Subject_ID);
E(E.Subject_ID==9,:) = [];

if ismember('HA_minus_RP_ROI9',E.Properties.VariableNames)
    E.HA_RP_Beta = E.HA_minus_RP_ROI9;
elseif all(ismember({'HA_Beta_ROI9','RP_Beta_ROI9'},E.Properties.VariableNames))
    E.HA_RP_Beta = E.HA_Beta_ROI9 - E.RP_Beta_ROI9;
else
    error(['ROI file must contain HA_minus_RP_ROI9 OR both ' ...
           'HA_Beta_ROI9 and RP_Beta_ROI9.']);
end

if numel(unique(E.Subject_ID)) ~= 57
    error('Expected N57 in ROI file; found N=%d.',numel(unique(E.Subject_ID)));
end

mu_beta = mean(E.HA_RP_Beta,'omitnan');
sd_beta = std(E.HA_RP_Beta,0,'omitnan');
E.Z_HA_RP_Beta = (E.HA_RP_Beta - mu_beta) ./ sd_beta;

P_GAIN = [0.40 0.50 0.60];
HIDDEN_MAG = [10 12.5 15];

nScen = numel(P_GAIN) * numel(HIDDEN_MAG);

Scenario = strings(nScen,1);
P_Gain = nan(nScen,1);
P_Loss = nan(nScen,1);
Hidden_Mag = nan(nScen,1);

LMM1_b = nan(nScen,1);
LMM1_SE = nan(nScen,1);
LMM1_CI_L = nan(nScen,1);
LMM1_CI_U = nan(nScen,1);
LMM1_F = nan(nScen,1);
LMM1_DF1 = nan(nScen,1);
LMM1_DF2 = nan(nScen,1);
LMM1_p = nan(nScen,1);

LMM2_b = nan(nScen,1);
LMM2_SE = nan(nScen,1);
LMM2_CI_L = nan(nScen,1);
LMM2_CI_U = nan(nScen,1);
LMM2_F = nan(nScen,1);
LMM2_DF1 = nan(nScen,1);
LMM2_DF2 = nan(nScen,1);
LMM2_p = nan(nScen,1);

LMM1_NTrials = nan(nScen,1);
LMM2_NTrials = nan(nScen,1);

row = 0;

for ip = 1:numel(P_GAIN)

    pGain = P_GAIN(ip);
    pLoss = 1 - pGain;

    for im = 1:numel(HIDDEN_MAG)

        row = row + 1;
        magAssumed = HIDDEN_MAG(im);

        fprintf('SCENARIO %d/%d: pGain=%.2f, pLoss=%.2f, hiddenMag=%.1f\n', ...
            row,nScen,pGain,pLoss,magAssumed);

        D = B;

        ambMask = strcmpi(D.Trial_Type,'Ambiguous');

        % Probability sensitivity.
        % Update ambiguous-trial probabilities.
        D.Gain_Prob(ambMask) = pGain;
        D.Loss_Prob(ambMask) = pLoss;

        % Magnitude sensitivity.
        % Only previously imputed 12.5 values are varied.
        gainHiddenMask = ambMask & abs(D.Gain_Mag-12.5)<1e-12;
        lossHiddenMask = ambMask & abs(D.Loss_Mag-12.5)<1e-12;

        D.Gain_Mag(gainHiddenMask) = magAssumed;
        D.Loss_Mag(lossHiddenMask) = magAssumed;

        D.SV_SENS = nan(height(D),1);

        for s = 1:height(Params)

            sid = Params.ID(s);
            idx = D.ID==sid;

            alpha  = Params.alpha(s);
            lambda = Params.lambda(s);
            gamma  = Params.gamma(s);

            pg = D.Gain_Prob(idx);
            pl = D.Loss_Prob(idx);
            gm = D.Gain_Mag(idx);
            lm = D.Loss_Mag(idx);

            wg = (pg.^gamma) ./ ...
                ((pg.^gamma + (1-pg).^gamma).^(1/gamma));

            wl = (pl.^gamma) ./ ...
                ((pl.^gamma + (1-pl).^gamma).^(1/gamma));

            D.SV_SENS(idx) = ...
                wg .* (gm.^alpha) - ...
                lambda .* wl .* (lm.^alpha);

        end

        D.Abs_SV_SENS = abs(D.SV_SENS);

        mu1 = mean(D.Abs_SV_SENS,'omitnan');
        sd1 = std(D.Abs_SV_SENS,0,'omitnan');
        D.Z_Abs_SV_LMM1 = (D.Abs_SV_SENS-mu1)./sd1;

        D38 = D(ismember(D.ID,EXTREME_IDS),:);

        grp38 = strings(height(D38),1);
        grp38(ismember(D38.ID,LOW_IDS)) = "Low CC";
        grp38(ismember(D38.ID,HIGH_IDS)) = "High CC";
        D38.CC_Group = categorical(grp38,{'Low CC','High CC'});
        D38.ID_cat = categorical(D38.ID);

        A38 = D38(strcmpi(D38.Trial_Type,'Ambiguous'),:);

        lmm1 = fitlme( ...
            A38, ...
            'Log_RT ~ Z_Abs_SV_LMM1 * CC_Group + (Z_Abs_SV_LMM1 | ID_cat)', ...
            'FitMethod','REML');

        A1 = anova(lmm1);

        s1 = extract_effect_stats( ...
            lmm1,A1,{'Z_Abs_SV_LMM1','CC_Group'});

        A57 = D(strcmpi(D.Trial_Type,'Ambiguous'),:);

        mu2 = mean(A57.Abs_SV_SENS,'omitnan');
        sd2 = std(A57.Abs_SV_SENS,0,'omitnan');
        A57.Z_Abs_SV_LMM2 = (A57.Abs_SV_SENS-mu2)./sd2;

        [matched,eidx] = ismember(A57.ID,E.Subject_ID);
        if any(~matched)
            error('Behavior-to-ROI merge failed.');
        end

        A57.Z_HA_RP_Beta = E.Z_HA_RP_Beta(eidx);

        A38_2 = A57(ismember(A57.ID,EXTREME_IDS),:);

        grp2 = strings(height(A38_2),1);
        grp2(ismember(A38_2.ID,LOW_IDS)) = "Low CC";
        grp2(ismember(A38_2.ID,HIGH_IDS)) = "High CC";

        A38_2.CC_Group = categorical(grp2,{'Low CC','High CC'});
        A38_2.ID_cat = categorical(A38_2.ID);

        lmm2 = fitlme( ...
            A38_2, ...
            ['Log_RT ~ Z_Abs_SV_LMM2 * Z_HA_RP_Beta * CC_Group ' ...
             '+ (1 | ID_cat)'], ...
            'FitMethod','REML');

        A2 = anova(lmm2);

        s2 = extract_effect_stats( ...
            lmm2,A2,{'Z_Abs_SV_LMM2','Z_HA_RP_Beta','CC_Group'});

        % Store results.

        Scenario(row) = sprintf('p%.2f_mag%.1f',pGain,magAssumed);
        P_Gain(row) = pGain;
        P_Loss(row) = pLoss;
        Hidden_Mag(row) = magAssumed;

        LMM1_b(row) = s1.b;
        LMM1_SE(row) = s1.SE;
        LMM1_CI_L(row) = s1.CI_L;
        LMM1_CI_U(row) = s1.CI_U;
        LMM1_F(row) = s1.F;
        LMM1_DF1(row) = s1.DF1;
        LMM1_DF2(row) = s1.DF2;
        LMM1_p(row) = s1.p;

        LMM2_b(row) = s2.b;
        LMM2_SE(row) = s2.SE;
        LMM2_CI_L(row) = s2.CI_L;
        LMM2_CI_U(row) = s2.CI_U;
        LMM2_F(row) = s2.F;
        LMM2_DF1(row) = s2.DF1;
        LMM2_DF2(row) = s2.DF2;
        LMM2_p(row) = s2.p;

        LMM1_NTrials(row) = height(A38);
        LMM2_NTrials(row) = height(A38_2);

        fprintf('LMM1 Ambiguous interaction: b=%.4f, F=%.4f, p=%.6g\n', ...
            s1.b,s1.F,s1.p);

        fprintf('LMM2 3-way interaction:     b=%.4f, F=%.4f, p=%.6g\n', ...
            s2.b,s2.F,s2.p);

    end
end

Results = table( ...
    Scenario,P_Gain,P_Loss,Hidden_Mag, ...
    LMM1_NTrials,LMM1_b,LMM1_SE,LMM1_CI_L,LMM1_CI_U, ...
    LMM1_F,LMM1_DF1,LMM1_DF2,LMM1_p, ...
    LMM2_NTrials,LMM2_b,LMM2_SE,LMM2_CI_L,LMM2_CI_U, ...
    LMM2_F,LMM2_DF1,LMM2_DF2,LMM2_p);

baseline = ...
    abs(Results.P_Gain-0.50)<1e-12 & ...
    abs(Results.Hidden_Mag-12.5)<1e-12;

if sum(baseline) ~= 1
    error('Could not identify the unique primary p=.50 / magnitude=12.5 scenario.');
end

b1_base = Results.LMM1_b(baseline);
b2_base = Results.LMM2_b(baseline);

Results.LMM1_SameDirection = sign(Results.LMM1_b)==sign(b1_base);
Results.LMM2_SameDirection = sign(Results.LMM2_b)==sign(b2_base);
Results.IsPrimaryAssumption = baseline;

disp(Results);

out_csv = fullfile(out_dir,'AMBIGUITY_PRIOR_SENSITIVITY_RESULTS.csv');
writetable(Results,out_csv);

save(fullfile(out_dir,'AMBIGUITY_PRIOR_SENSITIVITY_RESULTS.mat'), ...
    'Results','Params');

fprintf('Sensitivity analysis complete\n');
fprintf('Results saved to:\n%s\n',out_csv);

fprintf('\nPrimary assumption:\n');
disp(Results(baseline,:));

fprintf('\nRobustness summary:\n');
fprintf('LMM1 same-direction scenarios: %d / %d\n', ...
    sum(Results.LMM1_SameDirection),height(Results));
fprintf('LMM2 same-direction scenarios: %d / %d\n', ...
    sum(Results.LMM2_SameDirection),height(Results));
fprintf('LMM1 significant p<.05:        %d / %d\n', ...
    sum(Results.LMM1_p<.05),height(Results));
fprintf('LMM2 significant p<.05:        %d / %d\n', ...
    sum(Results.LMM2_p<.05),height(Results));

function nll = calc_NLL(p,lm,lp,gp,gm,ch)

    wg = (gp.^p(3)) ./ ...
        ((gp.^p(3) + (1-gp).^p(3)).^(1/p(3)));

    wl = (lp.^p(3)) ./ ...
        ((lp.^p(3) + (1-lp).^p(3)).^(1/p(3)));

    sv = wg .* (gm.^p(1)) - ...
         p(2) .* wl .* (lm.^p(1));

    pr = 1 ./ (1 + exp(-sv./p(4)));
    pr = max(min(pr,1-1e-10),1e-10);

    nll = -sum(ch.*log(pr) + (1-ch).*log(1-pr));

end

function S = extract_effect_stats(model,A,requiredTerms)

    names = string(model.Coefficients.Name);
    normNames = erase(names,[" ","(",")","{","}","'"]);

    mask = true(size(normNames));

    for k = 1:numel(requiredTerms)
        mask = mask & contains(normNames,requiredTerms{k});
    end

    % Match the requested interaction term.
    mask = mask & count(normNames,":")==numel(requiredTerms)-1;

    idxCoef = find(mask);

    if numel(idxCoef) ~= 1
        fprintf('\nCoefficient names:\n');
        disp(names);
        error('Could not uniquely identify target interaction coefficient.');
    end

    C = model.Coefficients;
    CI = coefCI(model);

    S.b = C.Estimate(idxCoef);

    % MATLAB compatibility.
    try
        S.SE = C.SE(idxCoef);
    catch
        V = model.CoefficientCovariance;
        S.SE = sqrt(V(idxCoef,idxCoef));
    end

    S.CI_L = CI(idxCoef,1);
    S.CI_U = CI(idxCoef,2);

    terms = get_anova_terms(A);
    normTerms = erase(terms,[" ","(",")","{","}","'"]);

    maskA = true(size(normTerms));

    for k = 1:numel(requiredTerms)
        maskA = maskA & contains(normTerms,requiredTerms{k});
    end

    maskA = maskA & count(normTerms,":")==numel(requiredTerms)-1;

    idxA = find(maskA);

    if numel(idxA) ~= 1
        fprintf('\nANOVA terms:\n');
        disp(terms);
        error('Could not uniquely identify target ANOVA interaction.');
    end

    S.F = double(A.FStat(idxA));
    S.DF1 = double(A.DF1(idxA));
    S.DF2 = double(A.DF2(idxA));
    S.p = double(A.pValue(idxA));

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
