clear;
clc;
close all;
%% EEGLAB path cleanup

bad_bcilab = ...
    'H:\eeglab13_5_4b\eeglab13_5_4b\plugins\Cleanline2.1\external\bcilab_partial';

if isfolder(bad_bcilab)

    bad_paths = strsplit(genpath(bad_bcilab),pathsep);
    current_paths = strsplit(path,pathsep);

    paths_to_remove = ...
        intersect(bad_paths,current_paths,'stable');

    for p = 1:numel(paths_to_remove)

        if ~isempty(paths_to_remove{p})
            rmpath(paths_to_remove{p});
        end

    end

end

clear Condition exp_condition;
rehash toolboxcache;
%% Paths

data_dir = 'H:\some cnt files copy\set files\ALL set files_MAIN\ultimate\final files';
save_dir = 'F:\my data\Figures\TF_MASTER_SIMPLE_FINAL';

if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% Analysis settings

epoch_window    = [-1500 4000];   % ms
baseline_window = [-200 0];       % ms
freq_range      = [13 30];        % Hz
roi_time        = [0 1000];       % ms
cycles          = [3 0.5];
MIN_TRIALS      = 3;

% A1/A2 are excluded from scalp analysis.
channels = {'FP1','FPz','FP2','F7','F3','Fz','F4','F8', ...
    'FC5','FC1','FC2','FC6','T7','C3','Cz','C4','T8', ...
    'CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8', ...
    'POz','O1','Oz','O2','AF7','AF3','AF4','AF8','F5','F1', ...
    'F2','F6','FC3','FCz','FC4','C5','C1','C2','C6','CP3', ...
    'CPz','CP4','P5','P1','P2','P6','PO5','PO3','PO4','PO6', ...
    'FT7','FT8','TP7','TP8','PO7','PO8'};

roi9 = {'F3','F4','Fz','FCz','P3','P4','Pz','CPz','POz'};
roi6 = {'F3','F4','Fz','P3','P4','Pz'};

nCh = numel(channels);
%% Event codes

rp_codes = arrayfun(@num2str,[16 18 19 34 36 37 43 46 48], ...
    'UniformOutput',false);

rn_codes = arrayfun(@num2str,[22 23 24 25 26 35 40 41 44], ...
    'UniformOutput',false);

la_codes = arrayfun(@num2str,[55:58 63:66], ...
    'UniformOutput',false);

ha_codes = arrayfun(@num2str,[59:62 67:72], ...
    'UniformOutput',false);

task_names = {'RP','RN','LA','HA'};
task_codes = {rp_codes,rn_codes,la_codes,ha_codes};
%% Subject files

files = dir(fullfile(data_dir,'*.set'));

if isempty(files)
    error('No .set files found in: %s',data_dir);
end

N = numel(files);
subject_id = NaN(N,1);

for s = 1:N
    tok = regexp(files(s).name,'^(\d+)_ref\.set$','tokens','once');

    if isempty(tok)
        error('Unexpected filename: %s',files(s).name);
    end

    subject_id(s) = str2double(tok{1}) - 100;
end

[subject_id,ord] = sort(subject_id);
files = files(ord);

fprintf('Subject files\n');
fprintf('N = %d\n',N);
disp(subject_id');

if numel(unique(subject_id)) ~= N
    error('Duplicate subject IDs found.');
end
%% Initialize

% grand stores RP, RN, LA, and HA in that order.
grand = cell(1,4);

times = [];
freqs = [];

subject_names = cell(N,1);
trial_counts = NaN(N,4);
task_ok = false(N,4);

log_id = [];
log_task = {};
log_channel = {};
log_problem = {};

% Use the first available channel location for each label.
chanlocs_selected = [];
have_chanloc = false(1,nCh);

fprintf('Starting TF analysis\n');
%% Time-frequency analysis

for s = 1:N

    sid = subject_id(s);
    filename = files(s).name;
    subject_names{s} = filename;

    fprintf('%d/%d | ID %d | %s\n',s,N,sid,filename);

    try
        EEG = pop_loadset('filename',filename,'filepath',data_dir);
        EEG = eeg_checkset(EEG);
    catch ME
        fprintf(2,'LOAD FAILED: %s\n',ME.message);

        for q = 1:4
            log_id(end+1,1) = sid;
            log_task{end+1,1} = task_names{q};
            log_channel{end+1,1} = 'ALL';
            log_problem{end+1,1} = ['LOAD FAILED: ' ME.message];
        end

        continue;
    end

    eeg_labels = {EEG.chanlocs.labels};

    if isempty(chanlocs_selected)
        chanlocs_selected = repmat(EEG.chanlocs(1),1,nCh);
        for ch = 1:nCh
            chanlocs_selected(ch).labels = channels{ch};
        end
    end

    for ch = 1:nCh
        idx_loc = find(strcmpi(eeg_labels,channels{ch}),1);

        if ~isempty(idx_loc) && ~have_chanloc(ch)
            chanlocs_selected(ch) = EEG.chanlocs(idx_loc);
            have_chanloc(ch) = true;
        end
    end

    present_event_types = cell(numel(EEG.event),1);

    for e = 1:numel(EEG.event)
        v = EEG.event(e).type;

        if isnumeric(v)
            present_event_types{e} = strtrim(num2str(v));
        elseif ischar(v)
            present_event_types{e} = strtrim(v);
        elseif isstring(v)
            present_event_types{e} = strtrim(char(v));
        else
            try
                present_event_types{e} = strtrim(char(v));
            catch
                present_event_types{e} = '';
            end
        end
    end

    present_event_types = unique(present_event_types);
%% Task loop

    for q = 1:4

        qname = task_names{q};
        requested = task_codes{q};

        % Use only event codes present in the dataset.
        present_requested = requested(ismember(requested,present_event_types));

        if isempty(present_requested)

            fprintf(2,'%s: no requested event markers found -> skipped.\n',qname);

            log_id(end+1,1) = sid;
            log_task{end+1,1} = qname;
            log_channel{end+1,1} = 'ALL';
            log_problem{end+1,1} = 'No requested event markers found';

            continue;
        end

        missing_requested = setdiff(requested,present_requested);

        if ~isempty(missing_requested)
            fprintf('%s: some event codes absent, continuing with present codes: ',qname);
            fprintf('%s ',present_requested{:});
            fprintf('\n');
        end

        try
            EEGq = pop_epoch(EEG,present_requested,epoch_window/1000);
            EEGq = eeg_checkset(EEGq);
        catch ME

            fprintf(2,'%s EPOCH FAILED: %s\n',qname,ME.message);

            log_id(end+1,1) = sid;
            log_task{end+1,1} = qname;
            log_channel{end+1,1} = 'ALL';
            log_problem{end+1,1} = ['EPOCH FAILED: ' ME.message];

            continue;
        end

        trial_counts(s,q) = EEGq.trials;

        fprintf('%s trials = %d\n',qname,EEGq.trials);

        if EEGq.trials < MIN_TRIALS

            fprintf(2,'%s: only %d trials -> skipped.\n',qname,EEGq.trials);

            log_id(end+1,1) = sid;
            log_task{end+1,1} = qname;
            log_channel{end+1,1} = 'ALL';
            log_problem{end+1,1} = sprintf('Too few trials: %d',EEGq.trials);

            continue;
        end

        try
            EEGq = pop_rmbase(EEGq,baseline_window);
        catch ME

            fprintf(2,'%s BASELINE FAILED: %s\n',qname,ME.message);

            log_id(end+1,1) = sid;
            log_task{end+1,1} = qname;
            log_channel{end+1,1} = 'ALL';
            log_problem{end+1,1} = ['BASELINE FAILED: ' ME.message];

            continue;
        end

        q_valid_channels = 0;

        for ch = 1:nCh

            cname = channels{ch};
            idx = find(strcmpi({EEGq.chanlocs.labels},cname));

            if isempty(idx)

                fprintf('  missing channel: %s\n',cname);

                log_id(end+1,1) = sid;
                log_task{end+1,1} = qname;
                log_channel{end+1,1} = cname;
                log_problem{end+1,1} = 'Missing channel';

                continue;
            end

            if numel(idx) ~= 1

                log_id(end+1,1) = sid;
                log_task{end+1,1} = qname;
                log_channel{end+1,1} = cname;
                log_problem{end+1,1} = sprintf('%d channel-label matches',numel(idx));

                continue;
            end

            try
                [ersp,~,~,t_out,f_out] = newtimef( ...
                    EEGq.data(idx,:,:), ...
                    EEGq.pnts, ...
                    epoch_window, ...
                    EEGq.srate, ...
                    cycles, ...
                    'baseline',baseline_window, ...
                    'freqs',freq_range, ...
                    'plotersp','off', ...
                    'plotitc','off', ...
                    'verbose','off');
            catch ME

                fprintf(2,'  newtimef failed: %s | %s\n',cname,ME.message);

                log_id(end+1,1) = sid;
                log_task{end+1,1} = qname;
                log_channel{end+1,1} = cname;
                log_problem{end+1,1} = ['newtimef failed: ' ME.message];

                continue;
            end

            if isempty(ersp) || any(~isfinite(ersp(:)))

                log_id(end+1,1) = sid;
                log_task{end+1,1} = qname;
                log_channel{end+1,1} = cname;
                log_problem{end+1,1} = 'Empty/NaN/Inf ERSP';

                continue;
            end

            % Set the common TF grid from the first successful channel.
            if isempty(times)

                times = t_out(:)';
                freqs = f_out(:)';

                nF = numel(freqs);
                nT = numel(times);

                for qq = 1:4
                    grand{qq} = NaN(nF,nT,nCh,N);
                end

                fprintf('\nTF grid established: %d freq x %d time bins\n',nF,nT);
                fprintf('Time %.3f to %.3f ms | Freq %.3f to %.3f Hz\n', ...
                    min(times),max(times),min(freqs),max(freqs));

            else

                if numel(t_out) ~= numel(times) || ...
                   numel(f_out) ~= numel(freqs) || ...
                   max(abs(double(t_out(:))-double(times(:)))) > 1e-8 || ...
                   max(abs(double(f_out(:))-double(freqs(:)))) > 1e-8

                    log_id(end+1,1) = sid;
                    log_task{end+1,1} = qname;
                    log_channel{end+1,1} = cname;
                    log_problem{end+1,1} = 'TF grid mismatch';

                    fprintf(2,'  TF grid mismatch: %s -> skipped.\n',cname);
                    continue;
                end
            end

            if ~isequal(size(ersp),[numel(freqs) numel(times)])

                log_id(end+1,1) = sid;
                log_task{end+1,1} = qname;
                log_channel{end+1,1} = cname;
                log_problem{end+1,1} = 'Unexpected ERSP dimensions';

                continue;
            end

            % Store only successful channel results.
            % Missing or failed channels remain NaN.
            grand{q}(:,:,ch,s) = ersp;

            q_valid_channels = q_valid_channels + 1;
        end

        if q_valid_channels > 0
            task_ok(s,q) = true;
        end

        fprintf('%s valid channels = %d/%d\n',qname,q_valid_channels,nCh);

    end
end
%% Check output

if isempty(times) || isempty(freqs)
    error('No successful TF result was produced.');
end

grand_ersp_RP = grand{1};
grand_ersp_RN = grand{2};
grand_ersp_LA = grand{3};
grand_ersp_HA = grand{4};

fprintf('TF dimensions\n');
disp(size(grand_ersp_RP));
%% Save logs

TrialCounts = table(subject_id, ...
    trial_counts(:,1),trial_counts(:,2),trial_counts(:,3),trial_counts(:,4), ...
    'VariableNames',{'Subject_ID','RP_Trials','RN_Trials','LA_Trials','HA_Trials'});

writetable(TrialCounts,fullfile(save_dir,'trial_counts.csv'));

TaskStatus = table(subject_id, ...
    task_ok(:,1),task_ok(:,2),task_ok(:,3),task_ok(:,4), ...
    'VariableNames',{'Subject_ID','RP_OK','RN_OK','LA_OK','HA_OK'});

writetable(TaskStatus,fullfile(save_dir,'task_status.csv'));

if isempty(log_id)

    ProblemLog = table();

else

    ProblemLog = table(log_id,log_task,log_channel,log_problem, ...
        'VariableNames',{'Subject_ID','Task','Channel','Problem'});

end

writetable(ProblemLog,fullfile(save_dir,'problem_log.csv'));
%% ROI indices

time_idx = find(times >= roi_time(1) & times <= roi_time(2));
freq_idx = find(freqs >= freq_range(1) & freqs <= freq_range(2));

fprintf('ROI bins\n');
fprintf('Time bins = %d | actual %.3f to %.3f ms\n', ...
    numel(time_idx),times(time_idx(1)),times(time_idx(end)));
fprintf('Freq bins = %d | actual %.3f to %.3f Hz\n', ...
    numel(freq_idx),freqs(freq_idx(1)),freqs(freq_idx(end)));
%% Channel-wise ROI values

nCols = 1 + nCh*4;
X = NaN(N,nCols);
X(:,1) = subject_id;

var_names = cell(1,nCols);
var_names{1} = 'Subject_ID';

col = 2;
for ch = 1:nCh
    var_names{col}   = [channels{ch} '_RP'];
    var_names{col+1} = [channels{ch} '_RN'];
    var_names{col+2} = [channels{ch} '_LA'];
    var_names{col+3} = [channels{ch} '_HA'];
    col = col + 4;
end

for s = 1:N

    col = 2;

    for ch = 1:nCh

        for q = 1:4

            x = grand{q}(freq_idx,time_idx,ch,s);
            finite_x = x(isfinite(x));

            if isempty(finite_x)
                val = NaN;
            else
                val = mean(finite_x);
            end

            X(s,col+q-1) = val;
        end

        col = col + 4;
    end
end

ROI_4Cond = array2table(X,'VariableNames',var_names);
ROI_4Cond = sortrows(ROI_4Cond,'Subject_ID');

roi_csv = fullfile(save_dir,'ROI_4Cond_pt_beta_FINAL.csv');
writetable(ROI_4Cond,roi_csv);
%% ROI summaries

RP9_cols = strcat(roi9,'_RP');
HA9_cols = strcat(roi9,'_HA');

RP6_cols = strcat(roi6,'_RP');
HA6_cols = strcat(roi6,'_HA');

RP9 = ROI_4Cond{:,RP9_cols};
HA9 = ROI_4Cond{:,HA9_cols};

RP6 = ROI_4Cond{:,RP6_cols};
HA6 = ROI_4Cond{:,HA6_cols};

complete_RP9 = all(isfinite(RP9),2);
complete_HA9 = all(isfinite(HA9),2);

complete_RP6 = all(isfinite(RP6),2);
complete_HA6 = all(isfinite(HA6),2);

RP_Beta_ROI9 = NaN(N,1);
HA_Beta_ROI9 = NaN(N,1);

RP_Beta_ROI6 = NaN(N,1);
HA_Beta_ROI6 = NaN(N,1);

RP_Beta_ROI9(complete_RP9) = mean(RP9(complete_RP9,:),2);
HA_Beta_ROI9(complete_HA9) = mean(HA9(complete_HA9,:),2);

RP_Beta_ROI6(complete_RP6) = mean(RP6(complete_RP6,:),2);
HA_Beta_ROI6(complete_HA6) = mean(HA6(complete_HA6,:),2);

HA_minus_RP_ROI9 = HA_Beta_ROI9 - RP_Beta_ROI9;
HA_minus_RP_ROI6 = HA_Beta_ROI6 - RP_Beta_ROI6;

ROI9_Complete = isfinite(RP_Beta_ROI9) & isfinite(HA_Beta_ROI9);
ROI6_Complete = isfinite(RP_Beta_ROI6) & isfinite(HA_Beta_ROI6);

ROI_Summary = table(subject_id, ...
    RP_Beta_ROI9,HA_Beta_ROI9,HA_minus_RP_ROI9,ROI9_Complete, ...
    RP_Beta_ROI6,HA_Beta_ROI6,HA_minus_RP_ROI6,ROI6_Complete, ...
    'VariableNames',{'Subject_ID', ...
    'RP_Beta_ROI9','HA_Beta_ROI9','HA_minus_RP_ROI9','ROI9_Complete', ...
    'RP_Beta_ROI6','HA_Beta_ROI6','HA_minus_RP_ROI6','ROI6_Complete'});

ROI_Summary = sortrows(ROI_Summary,'Subject_ID');

writetable(ROI_Summary,fullfile(save_dir,'ROI9_ROI6_summary_FINAL.csv'));
%% Paired ROI9 test

mask = ROI9_Complete;

fprintf('ROI9 completeness\n');
fprintf('Complete ROI9 RP+HA = %d/%d\n',sum(mask),N);

if any(~mask)
    fprintf('Incomplete ROI9 IDs:\n');
    disp(subject_id(~mask)');
end

if sum(mask) >= 2

    d = HA_minus_RP_ROI9(mask);

    [~,p,ci,st] = ttest(d,0);

    dz = mean(d) / std(d,0);

    fprintf('ROI9 statistics\n');

    fprintf('N = %d\n',sum(mask));
    fprintf('RP mean = %.9f | SD = %.9f\n', ...
        mean(RP_Beta_ROI9(mask)),std(RP_Beta_ROI9(mask),0));

    fprintf('HA mean = %.9f | SD = %.9f\n', ...
        mean(HA_Beta_ROI9(mask)),std(HA_Beta_ROI9(mask),0));

    fprintf('HA-RP = %.9f | SD = %.9f\n', ...
        mean(d),std(d,0));

    fprintf('t(%d) = %.9f | p = %.12g\n', ...
        st.df,st.tstat,p);

    fprintf('95%% CI = [%.9f, %.9f]\n',ci(1),ci(2));
    fprintf('Cohen dz = %.9f\n',dz);
end
%% Save MAT file

mat_file = fullfile(save_dir,'TF_OUTPUT_FINAL.mat');

save(mat_file, ...
    'grand_ersp_RP','grand_ersp_RN','grand_ersp_LA','grand_ersp_HA', ...
    'times','freqs','channels','chanlocs_selected','have_chanloc', ...
    'subject_names','subject_id','trial_counts','task_ok', ...
    'roi9','roi6','epoch_window','baseline_window','freq_range','roi_time','cycles', ...
    '-v7.3');
%% Summary

fprintf('Done\n');

fprintf('RP usable task data: %d/%d\n',sum(task_ok(:,1)),N);
fprintf('RN usable task data: %d/%d\n',sum(task_ok(:,2)),N);
fprintf('LA usable task data: %d/%d\n',sum(task_ok(:,3)),N);
fprintf('HA usable task data: %d/%d\n',sum(task_ok(:,4)),N);
fprintf('Strict ROI9 complete: %d/%d\n',sum(ROI9_Complete),N);

fprintf('\nMAT:\n%s\n',mat_file);
fprintf('\nROI CSV:\n%s\n',roi_csv);
fprintf('\nProblem log:\n%s\n',fullfile(save_dir,'problem_log.csv'));

fprintf('\nMissing/failed values stayed NaN; none were replaced by zero.\n');
