clear;
clc;
close all;

data_dir = ...
    'H:\some cnt files copy\set files\ALL set files_MAIN\ultimate\final files';

primary_roi_csv = ...
    'F:\my data\ROI_4Cond_pt_beta_all2.csv';

save_dir = ...
    'F:\my data\Figures\HA_RP_RELIABILITY_EXACT_FINAL';

if ~isfolder(data_dir)
    error('DATA DIRECTORY NOT FOUND:\n%s',data_dir);
end

if ~isfile(primary_roi_csv)
    error(['PRIMARY ROI CSV NOT FOUND:\n%s\n\n' ...
           'This should be the ROI_4Cond_pt_beta_FINAL.csv exported by ' ...
           'the MASTER TF pipeline.'], ...
          primary_roi_csv);
end

if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

required_functions = { ...
    'pop_loadset', ...
    'eeg_checkset', ...
    'pop_epoch', ...
    'pop_rmbase', ...
    'pop_select', ...
    'newtimef'};

for k = 1:numel(required_functions)
    if exist(required_functions{k},'file') ~= 2
        error('Required EEGLAB function is not on MATLAB path: %s', ...
            required_functions{k});
    end
end

fprintf('EEGLAB functions\n');
fprintf('All required EEGLAB functions are available.\n');

epoch_window = [-1500 4000];        % ms
baseline_window = [-200 0];         % ms
freq_range = [13 30];               % Hz
analysis_window = [0 1000];         % ms
wavelet_cycles = [3 0.5];

roi9 = { ...
    'F3','F4','Fz','FCz', ...
    'P3','P4','Pz','CPz','POz'};

N_SPLITS = 50;

RNG_SEED = 20260916;
MIN_TRIALS_FOR_RELIABILITY = 4;

rng(RNG_SEED,'twister');

subj_pos_markers = [16 18 19 34 36 37 43 46 48];
high_amb_markers = [59:62 67:72];

ev_rp = arrayfun(@num2str,subj_pos_markers,'UniformOutput',false);
ev_ha = arrayfun(@num2str,high_amb_markers,'UniformOutput',false);

EXPECTED_IDS = [ ...
     1  2  3  4  6  8 11 12 13 14 ...
    15 16 17 18 19 20 21 22 23 25 ...
    26 27 28 29 31 32 33 36 37 38 ...
    39 40 41 42 43 45 46 47 48 49 ...
    50 52 53 54 55 56 57 58 60 61 ...
    63 64 65 66 67 68 69];

LOW_IDS = [ ...
     1  2 11 12 13 15 27 28 38 40 ...
    45 47 48 50 52 53 60 66 67];

HIGH_IDS = [ ...
     6 17 21 22 25 29 31 32 33 36 ...
    37 39 41 46 49 54 55 58 61];

EXTREME38_IDS = sort([LOW_IDS HIGH_IDS]);

if numel(EXPECTED_IDS) ~= 57 || numel(unique(EXPECTED_IDS)) ~= 57
    error('EXPECTED_IDS is not exactly 57 unique IDs.');
end

if numel(EXTREME38_IDS) ~= 38 || numel(unique(EXTREME38_IDS)) ~= 38
    error('EXTREME38_IDS is not exactly 38 unique IDs.');
end

PrimaryTable = readtable(primary_roi_csv);

required_RP_cols = strcat(roi9,'_RP');
required_HA_cols = strcat(roi9,'_HA');
required_cols = [{'Subject_ID'} required_RP_cols required_HA_cols];

missing_cols = setdiff(required_cols,PrimaryTable.Properties.VariableNames);

if ~isempty(missing_cols)
    fprintf('\nMissing required columns from primary ROI CSV:\n');
    disp(missing_cols');
    error('Primary ROI CSV does not contain all required ROI9 columns.');
end

PrimaryTable = sortrows(PrimaryTable,'Subject_ID');

if height(PrimaryTable) ~= 57
    error('Primary ROI CSV must contain 57 rows; found %d.',height(PrimaryTable));
end

if ~isequal(PrimaryTable.Subject_ID(:)',sort(EXPECTED_IDS))
    error('Primary ROI CSV does not contain the exact expected N57 IDs.');
end

Primary_RP_ROI9 = mean(PrimaryTable{:,required_RP_cols},2);
Primary_HA_ROI9 = mean(PrimaryTable{:,required_HA_cols},2);
Primary_DIFF_ROI9 = Primary_HA_ROI9 - Primary_RP_ROI9;

fprintf('Primary ROI data loaded\n');
fprintf('N = %d\n',height(PrimaryTable));
fprintf('RP ROI9 mean    = %.9f dB\n',mean(Primary_RP_ROI9));
fprintf('HA ROI9 mean    = %.9f dB\n',mean(Primary_HA_ROI9));
fprintf('HA-RP ROI9 mean = %.9f dB\n',mean(Primary_DIFF_ROI9));

subject_files = dir(fullfile(data_dir,'*.set'));

if isempty(subject_files)
    error('No .set files found in:\n%s',data_dir);
end

n_files = numel(subject_files);
raw_names = string({subject_files.name})';

EEG_ID = NaN(n_files,1);
Subject_ID = NaN(n_files,1);

for s = 1:n_files
    [EEG_ID(s),Subject_ID(s)] = parse_subject_filename(raw_names(s));
end

if numel(unique(Subject_ID)) ~= n_files
    error('Duplicate Subject_ID values detected.');
end

[Subject_ID,sort_idx] = sort(Subject_ID);
EEG_ID = EEG_ID(sort_idx);
subject_files = subject_files(sort_idx);
raw_names = raw_names(sort_idx);

if n_files ~= 57
    error('Expected exactly 57 final .set files; found %d.',n_files);
end

if ~isequal(Subject_ID(:)',sort(EXPECTED_IDS))
    error('Final .set subject pool does not match expected N57.');
end

if ~isequal(Subject_ID(:),PrimaryTable.Subject_ID(:))
    error('Subject order in .set files does not match primary ROI CSV.');
end

N = n_files;
extreme38_mask = ismember(Subject_ID,EXTREME38_IDS);

if sum(extreme38_mask) ~= 38
    error('Could not recover the exact N38 extreme-groups sample.');
end

nRP = NaN(N,1);
nHA = NaN(N,1);

Full_RP_Reconstructed = NaN(N,1);
Full_HA_Reconstructed = NaN(N,1);

RP_half1 = NaN(N,N_SPLITS);
RP_half2 = NaN(N,N_SPLITS);
HA_half1 = NaN(N,N_SPLITS);
HA_half2 = NaN(N,N_SPLITS);
DIFF_half1 = NaN(N,N_SPLITS);
DIFF_half2 = NaN(N,N_SPLITS);

fprintf('Split-half reliability analysis\n');
fprintf('N subjects = %d\n',N);
fprintf('N random splits = %d\n',N_SPLITS);

for s = 1:N

    filename = subject_files(s).name;
    this_id = Subject_ID(s);

    fprintf('Subject %d/%d | ID %d | %s\n',s,N,this_id,filename);

    EEG = pop_loadset('filename',filename,'filepath',data_dir);
    EEG = eeg_checkset(EEG);

    labels = string({EEG.chanlocs.labels});

    for r = 1:numel(roi9)
        idx = find(strcmpi(labels,roi9{r}));
        if numel(idx) ~= 1
            error('Subject %d: ROI9 channel lookup failed for %s.', ...
                this_id,roi9{r});
        end
    end

    EEG_RP = pop_epoch(EEG,ev_rp,epoch_window/1000);
    EEG_HA = pop_epoch(EEG,ev_ha,epoch_window/1000);

    EEG_RP = pop_rmbase(EEG_RP,baseline_window);
    EEG_HA = pop_rmbase(EEG_HA,baseline_window);

    nRP(s) = EEG_RP.trials;
    nHA(s) = EEG_HA.trials;

    fprintf('Trials: RP=%d | HA=%d\n',nRP(s),nHA(s));

    if nRP(s) < MIN_TRIALS_FOR_RELIABILITY || ...
       nHA(s) < MIN_TRIALS_FOR_RELIABILITY
        error(['Subject %d has too few trials for split-half reliability. ' ...
               'RP=%d, HA=%d.'],this_id,nRP(s),nHA(s));
    end

    % Reconstruct full-condition ROI9 values.
    Full_RP_Reconstructed(s) = compute_roi9_beta_exact( ...
        EEG_RP,roi9,epoch_window,baseline_window, ...
        freq_range,wavelet_cycles,analysis_window);

    Full_HA_Reconstructed(s) = compute_roi9_beta_exact( ...
        EEG_HA,roi9,epoch_window,baseline_window, ...
        freq_range,wavelet_cycles,analysis_window);

    dRP = abs(Full_RP_Reconstructed(s)-Primary_RP_ROI9(s));
    dHA = abs(Full_HA_Reconstructed(s)-Primary_HA_ROI9(s));

    fprintf('Full-condition check: RP diff=%.3g | HA diff=%.3g\n',dRP,dHA);

    for split_i = 1:N_SPLITS

        [rp_idx1,rp_idx2] = balanced_random_halves(nRP(s));
        [ha_idx1,ha_idx2] = balanced_random_halves(nHA(s));

        EEG_RP_1 = pop_select(EEG_RP,'trial',rp_idx1);
        EEG_RP_2 = pop_select(EEG_RP,'trial',rp_idx2);
        EEG_HA_1 = pop_select(EEG_HA,'trial',ha_idx1);
        EEG_HA_2 = pop_select(EEG_HA,'trial',ha_idx2);

        RP_half1(s,split_i) = compute_roi9_beta_exact( ...
            EEG_RP_1,roi9,epoch_window,baseline_window, ...
            freq_range,wavelet_cycles,analysis_window);

        RP_half2(s,split_i) = compute_roi9_beta_exact( ...
            EEG_RP_2,roi9,epoch_window,baseline_window, ...
            freq_range,wavelet_cycles,analysis_window);

        HA_half1(s,split_i) = compute_roi9_beta_exact( ...
            EEG_HA_1,roi9,epoch_window,baseline_window, ...
            freq_range,wavelet_cycles,analysis_window);

        HA_half2(s,split_i) = compute_roi9_beta_exact( ...
            EEG_HA_2,roi9,epoch_window,baseline_window, ...
            freq_range,wavelet_cycles,analysis_window);

        DIFF_half1(s,split_i) = HA_half1(s,split_i)-RP_half1(s,split_i);
        DIFF_half2(s,split_i) = HA_half2(s,split_i)-RP_half2(s,split_i);

        if mod(split_i,10)==0 || split_i==N_SPLITS
            fprintf('  Completed split %d/%d\n',split_i,N_SPLITS);
        end
    end
end

max_rp_discrepancy = max(abs(Full_RP_Reconstructed-Primary_RP_ROI9));
max_ha_discrepancy = max(abs(Full_HA_Reconstructed-Primary_HA_ROI9));
max_diff_discrepancy = max(abs( ...
    (Full_HA_Reconstructed-Full_RP_Reconstructed)-Primary_DIFF_ROI9));

fprintf('Reconstruction check\n');
fprintf('Max RP discrepancy    = %.15g dB\n',max_rp_discrepancy);
fprintf('Max HA discrepancy    = %.15g dB\n',max_ha_discrepancy);
fprintf('Max HA-RP discrepancy = %.15g dB\n',max_diff_discrepancy);

RECON_TOL = 1e-10;

if max_rp_discrepancy > RECON_TOL || max_ha_discrepancy > RECON_TOL
    error(['EXACT RECONSTRUCTION FAILED. Do not interpret reliability. ' ...
           'The full-condition NEWTIMEF values do not reproduce the ' ...
           'primary ROI CSV.']);
end

fprintf('Full-condition reconstruction matches the primary pipeline.\n');

[rDiff57_raw,rDiff57_SB] = repeated_reliability( ...
    DIFF_half1,DIFF_half2,true(N,1));

[rHA57_raw,rHA57_SB] = repeated_reliability( ...
    HA_half1,HA_half2,true(N,1));

[rRP57_raw,rRP57_SB] = repeated_reliability( ...
    RP_half1,RP_half2,true(N,1));

[rDiff38_raw,rDiff38_SB] = repeated_reliability( ...
    DIFF_half1,DIFF_half2,extreme38_mask);

[rHA38_raw,rHA38_SB] = repeated_reliability( ...
    HA_half1,HA_half2,extreme38_mask);

[rRP38_raw,rRP38_SB] = repeated_reliability( ...
    RP_half1,RP_half2,extreme38_mask);

Summary57 = summarize_reliability( ...
    'N57',rDiff57_raw,rDiff57_SB, ...
    rHA57_raw,rHA57_SB,rRP57_raw,rRP57_SB);

Summary38 = summarize_reliability( ...
    'N38',rDiff38_raw,rDiff38_SB, ...
    rHA38_raw,rHA38_SB,rRP38_raw,rRP38_SB);

ReliabilitySummary = [Summary57; Summary38];

fprintf('Reliability results\n');
disp(ReliabilitySummary);

row57 = strcmp(ReliabilitySummary.Sample,'N57') & ...
        strcmp(ReliabilitySummary.Measure,'HA-RP difference');

row38 = strcmp(ReliabilitySummary.Sample,'N38') & ...
        strcmp(ReliabilitySummary.Measure,'HA-RP difference');

fprintf('\nN57:\n');
fprintf('HA-RP median Spearman-Brown r = %.3f\n', ...
    ReliabilitySummary.Median_SB_r(row57));
fprintf('2.5th-97.5th percentile = [%.3f, %.3f]\n', ...
    ReliabilitySummary.P2_5_SB(row57), ...
    ReliabilitySummary.P97_5_SB(row57));

fprintf('\nN38:\n');
fprintf('HA-RP median Spearman-Brown r = %.3f\n', ...
    ReliabilitySummary.Median_SB_r(row38));
fprintf('2.5th-97.5th percentile = [%.3f, %.3f]\n', ...
    ReliabilitySummary.P2_5_SB(row38), ...
    ReliabilitySummary.P97_5_SB(row38));

summary_csv = fullfile(save_dir,'HA_RP_reliability_summary_EXACT.csv');
writetable(ReliabilitySummary,summary_csv);

SplitResults = table( ...
    (1:N_SPLITS)', ...
    rDiff57_raw,rDiff57_SB, ...
    rHA57_raw,rHA57_SB, ...
    rRP57_raw,rRP57_SB, ...
    rDiff38_raw,rDiff38_SB, ...
    rHA38_raw,rHA38_SB, ...
    rRP38_raw,rRP38_SB, ...
    'VariableNames',{ ...
    'Split', ...
    'N57_DIFF_raw_r','N57_DIFF_SB_r', ...
    'N57_HA_raw_r','N57_HA_SB_r', ...
    'N57_RP_raw_r','N57_RP_SB_r', ...
    'N38_DIFF_raw_r','N38_DIFF_SB_r', ...
    'N38_HA_raw_r','N38_HA_SB_r', ...
    'N38_RP_raw_r','N38_RP_SB_r'});

split_csv = fullfile(save_dir,'HA_RP_reliability_all_splits_EXACT.csv');
writetable(SplitResults,split_csv);

SubjectAudit = table( ...
    Subject_ID,nRP,nHA, ...
    Primary_RP_ROI9,Full_RP_Reconstructed, ...
    Primary_HA_ROI9,Full_HA_Reconstructed, ...
    Primary_DIFF_ROI9, ...
    Full_HA_Reconstructed-Full_RP_Reconstructed, ...
    extreme38_mask, ...
    'VariableNames',{ ...
    'Subject_ID','RP_Trials','HA_Trials', ...
    'Primary_RP_ROI9','Reconstructed_RP_ROI9', ...
    'Primary_HA_ROI9','Reconstructed_HA_ROI9', ...
    'Primary_HAminusRP_ROI9','Reconstructed_HAminusRP_ROI9', ...
    'In_Extreme38'});

audit_csv = fullfile(save_dir,'HA_RP_full_condition_reconstruction_EXACT.csv');
writetable(SubjectAudit,audit_csv);

mat_output = fullfile(save_dir,'HA_RP_reliability_EXACT.mat');

save(mat_output, ...
    'ReliabilitySummary','SplitResults','SubjectAudit', ...
    'RP_half1','RP_half2','HA_half1','HA_half2', ...
    'DIFF_half1','DIFF_half2','Subject_ID','EXTREME38_IDS', ...
    'N_SPLITS','RNG_SEED', ...
    'max_rp_discrepancy','max_ha_discrepancy','max_diff_discrepancy', ...
    '-v7.3');

fprintf('\nOutputs saved to:\n%s\n',save_dir);
fprintf('\nDone.\n');

function beta_roi9 = compute_roi9_beta_exact( ...
    EEGc,roi9,epoch_window,baseline_window, ...
    freq_range,wavelet_cycles,analysis_window)

    labels = string({EEGc.chanlocs.labels});
    channel_beta = NaN(numel(roi9),1);

    ref_times = [];
    ref_freqs = [];

    for r = 1:numel(roi9)

        target_chan = roi9{r};
        idx = find(strcmpi(labels,target_chan));

        if numel(idx) ~= 1
            error('ROI9 channel lookup failed for %s.',target_chan);
        end

        [this_ersp,~,~,this_times,this_freqs] = ...
            newtimef( ...
                EEGc.data(idx,:,:), ...
                EEGc.pnts, ...
                epoch_window, ...
                EEGc.srate, ...
                wavelet_cycles, ...
                'baseline',baseline_window, ...
                'freqs',freq_range, ...
                'plotersp','off', ...
                'plotitc','off', ...
                'verbose','off');

        if isempty(ref_times)
            ref_times = this_times(:)';
            ref_freqs = this_freqs(:)';
        else
            assert_same_numeric_vector(ref_times,this_times, ...
                sprintf('ROI9 time grid at %s',target_chan));
            assert_same_numeric_vector(ref_freqs,this_freqs, ...
                sprintf('ROI9 frequency grid at %s',target_chan));
        end

        time_idx = find(this_times >= analysis_window(1) & ...
                        this_times <= analysis_window(2));

        freq_idx = find(this_freqs >= freq_range(1) & ...
                        this_freqs <= freq_range(2));

        if isempty(time_idx)
            error('No NEWTIMEF time bins lie in the analysis window.');
        end

        if isempty(freq_idx)
            error('No NEWTIMEF frequency bins lie in the beta range.');
        end

        x = this_ersp(freq_idx,time_idx);
        channel_beta(r) = mean(x(:));

    end

    if any(~isfinite(channel_beta))
        error('Non-finite ROI9 channel beta estimate.');
    end

    beta_roi9 = mean(channel_beta);

end

function [idx1,idx2] = balanced_random_halves(nTrials)

    if nTrials < 4
        error('At least 4 trials are required for split-half reliability.');
    end

    nHalf = floor(nTrials/2);
    p = randperm(nTrials);

    % Use equal-sized halves.
    % If the trial count is odd, one trial is omitted.
    idx1 = p(1:nHalf);
    idx2 = p(nHalf+1:2*nHalf);

end

function [rRaw,rSB] = repeated_reliability(X1,X2,subject_mask)

    nSplits = size(X1,2);
    rRaw = NaN(nSplits,1);
    rSB = NaN(nSplits,1);

    subject_mask = logical(subject_mask(:));

    for k = 1:nSplits

        a = X1(subject_mask,k);
        b = X2(subject_mask,k);

        good = isfinite(a) & isfinite(b);

        if sum(good) < 3
            continue;
        end

        r = corr(a(good),b(good), ...
            'Type','Pearson','Rows','complete');

        rRaw(k) = r;
        rSB(k) = spearman_brown(r);

    end

end

function T = summarize_reliability( ...
    sample_name,rDiffRaw,rDiffSB,rHARaw,rHASB,rRPRaw,rRPSB)

    Measure = {'HA-RP difference';'HA';'RP'};
    Sample = repmat({sample_name},3,1);

    Median_Raw_r = [ ...
        median(rDiffRaw(isfinite(rDiffRaw))); ...
        median(rHARaw(isfinite(rHARaw))); ...
        median(rRPRaw(isfinite(rRPRaw)))];

    Median_SB_r = [ ...
        median(rDiffSB(isfinite(rDiffSB))); ...
        median(rHASB(isfinite(rHASB))); ...
        median(rRPSB(isfinite(rRPSB)))];

    ciDiff = prctile(rDiffSB(isfinite(rDiffSB)),[2.5 97.5]);
    ciHA = prctile(rHASB(isfinite(rHASB)),[2.5 97.5]);
    ciRP = prctile(rRPSB(isfinite(rRPSB)),[2.5 97.5]);

    P2_5_SB = [ciDiff(1);ciHA(1);ciRP(1)];
    P97_5_SB = [ciDiff(2);ciHA(2);ciRP(2)];

    T = table(Sample,Measure,Median_Raw_r,Median_SB_r,P2_5_SB,P97_5_SB);

end

function rSB = spearman_brown(r)

    if ~isfinite(r) || abs(1+r) < eps
        rSB = NaN;
    else
        rSB = (2*r)/(1+r);
    end

end

function [eeg_id,subject_id] = parse_subject_filename(filename)

    filename = strtrim(char(filename));

    token = regexp(filename,'^(\d+)_ref\.set$','tokens','once');

    if isempty(token)
        error(['Filename does not match required pattern ' ...
               '"###_ref.set": %s'],filename);
    end

    eeg_id = str2double(token{1});
    subject_id = eeg_id - 100;

    if ~isfinite(eeg_id) || ~isfinite(subject_id)
        error('Invalid numeric ID parsed from %s.',filename);
    end

end

function assert_same_numeric_vector(a,b,label)

    a = double(a(:));
    b = double(b(:));

    if numel(a) ~= numel(b)
        error('%s: vector lengths differ.',label);
    end

    tol = 1e-10;

    if max(abs(a-b)) > tol
        error('%s: maximum discrepancy %.15g exceeds tolerance %.3g.', ...
            label,max(abs(a-b)),tol);
    end

end
