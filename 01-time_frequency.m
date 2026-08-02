% =========================================================================
% TFA: Risky (EV+) vs. Risky (EV-) vs. Low Ambiguity vs. High Ambiguity
% Pure Stimulus-Locked Valuation - No Trial Equalization
% =========================================================================
clear; clc;


%% 1.
data_dir = 'H:\some cnt files copy\set files\ALL set files_MAIN\ultimate';
save_dir = 'F:\my data\Figures\';
channels_of_interest = {'FP1', 'FPz', 'FP2', 'F7', 'F3', 'Fz', 'F4', 'F8', ...
    'FC5', 'FC1', 'FC2', 'FC6', 'A1', 'T7', 'C3', 'Cz', 'C4', 'T8', ...
    'A2', 'CP5', 'CP1', 'CP2', 'CP6', 'P7', 'P3', 'Pz', 'P4', 'P8', ...
    'POz', 'O1', 'Oz', 'O2', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', ...
    'F2', 'F6', 'FC3', 'FCz', 'FC4', 'C5', 'C1', 'C2', 'C6', 'CP3', ...
    'CPz', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', ...
    'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8'};
% epoch and baseline
epoch_window    = [-1500 4000];
baseline_window = [-200 0];
freq_range      = [13 30];

%% 2. markers interpretation
subj_pos_markers = [16, 18, 19, 34, 36, 37, 43, 46, 48];
subj_neg_markers = [22, 23, 24, 25, 26, 35, 40, 41, 44];
low_amb_markers  = [55:58, 63:66]; 
high_amb_markers = [59:62, 67:72];

% 
ev_rp = strtrim(cellstr(num2str(subj_pos_markers')));
ev_rn = strtrim(cellstr(num2str(subj_neg_markers')));
ev_la = strtrim(cellstr(num2str(low_amb_markers')));
ev_ha = strtrim(cellstr(num2str(high_amb_markers')));

grand_ersp_RP = [];
grand_ersp_RN = [];
grand_ersp_LA = [];
grand_ersp_HA = [];

valid_subjects = 0;
times = []; freqs = [];

subject_files = dir(fullfile(data_dir, '*.set'));

%% 3. process
for s = 1:length(subject_files)
    filename = subject_files(s).name;
    fprintf('\n---> Processing Subject: %s\n', filename);
    
    EEG = pop_loadset('filename', filename, 'filepath', data_dir);
    
    try
        % epoching
        EEG_RP = pop_epoch(EEG, ev_rp, epoch_window/1000);
        EEG_RP = pop_rmbase(EEG_RP, baseline_window);
        
        EEG_RN = pop_epoch(EEG, ev_rn, epoch_window/1000);
        EEG_RN = pop_rmbase(EEG_RN, baseline_window);
        
        EEG_LA = pop_epoch(EEG, ev_la,  epoch_window/1000);
        EEG_LA = pop_rmbase(EEG_LA, baseline_window);
        
        EEG_HA = pop_epoch(EEG, ev_ha,  epoch_window/1000);
        EEG_HA = pop_rmbase(EEG_HA, baseline_window);

        % assessing trial numbers
        if EEG_RP.trials > 2 && EEG_RN.trials > 2 && EEG_LA.trials > 2 && EEG_HA.trials > 2
            
            temp_rp = []; temp_rn = []; temp_la = []; temp_ha = [];
            subj_success = true;
            
            for ch = 1:length(channels_of_interest)
                target_chan = channels_of_interest{ch};
                chan_idx = find(strcmpi({EEG.chanlocs.labels}, target_chan));
                
                if ~isempty(chan_idx)
                    try
                        [ersp_rp, ~, ~, t_out, f_out] = newtimef(EEG_RP.data(chan_idx,:,:), ...
                            EEG_RP.pnts, epoch_window, EEG_RP.srate, [3 0.5], 'baseline', baseline_window, ...
                            'freqs', freq_range, 'plotersp', 'off', 'plotitc', 'off');
                            
                        [ersp_rn] = newtimef(EEG_RN.data(chan_idx,:,:), ...
                            EEG_RN.pnts, epoch_window, EEG_RN.srate, [3 0.5], 'baseline', baseline_window, ...
                            'freqs', freq_range, 'plotersp', 'off', 'plotitc', 'off');
                            
                        [ersp_la] = newtimef(EEG_LA.data(chan_idx,:,:), ...
                            EEG_LA.pnts, epoch_window, EEG_LA.srate, [3 0.5], 'baseline', baseline_window, ...
                            'freqs', freq_range, 'plotersp', 'off', 'plotitc', 'off');
                            
                        [ersp_ha] = newtimef(EEG_HA.data(chan_idx,:,:), ...
                            EEG_HA.pnts, epoch_window, EEG_HA.srate, [3 0.5], 'baseline', baseline_window, ...
                            'freqs', freq_range, 'plotersp', 'off', 'plotitc', 'off');
                            
                        temp_rp(:,:,ch) = ersp_rp;
                        temp_rn(:,:,ch) = ersp_rn;
                        temp_la(:,:,ch) = ersp_la;
                        temp_ha(:,:,ch) = ersp_ha;
                        
                        times = t_out; freqs = f_out; 
                    catch
                        subj_success = false; 
                        break;
                    end
                end
            end
            
            if subj_success && ~isempty(temp_rp)
                valid_subjects = valid_subjects + 1;
                grand_ersp_RP(:,:,:,valid_subjects) = temp_rp;
                grand_ersp_RN(:,:,:,valid_subjects) = temp_rn;
                grand_ersp_LA(:,:,:,valid_subjects) = temp_la;
                grand_ersp_HA(:,:,:,valid_subjects) = temp_ha;
            end
        end
    catch
        continue; 
    end
end

%% =========================================================================
% Step 2: ROI Data Extraction for 4 Conditions
% =========================================================================
roi_time = [0 1000];            
roi_freq = [13 30];              
save_excel_path = 'F:\my data\Figures\ROI_4Cond_pt_beta_all2_1.xlsx';

time_idx = find(times >= roi_time(1) & times <= roi_time(2));
freq_idx = find(freqs >= roi_freq(1) & freqs <= roi_freq(2));

if isempty(time_idx) || isempty(freq_idx)
    error('Error! the frequency or time is out of range');
end

num_subj = size(grand_ersp_RP, 4);
num_chans = length(channels_of_interest);

extracted_data = zeros(num_subj, 1 + (num_chans * 4));
extracted_data(:, 1) = 1:num_subj;

var_names = {'Subject_ID'};
for ch = 1:num_chans
    chan_name = channels_of_interest{ch};
    var_names{end+1} = sprintf('%s_RP', chan_name); 
    var_names{end+1} = sprintf('%s_RN', chan_name); 
    var_names{end+1} = sprintf('%s_LA', chan_name); % Low Ambiguity
    var_names{end+1} = sprintf('%s_HA', chan_name); % High Ambiguity
end

for s = 1:num_subj
    col_idx = 2;
    for ch = 1:num_chans
        subj_rp = squeeze(grand_ersp_RP(freq_idx, time_idx, ch, s));
        subj_rn = squeeze(grand_ersp_RN(freq_idx, time_idx, ch, s));
        subj_la = squeeze(grand_ersp_LA(freq_idx, time_idx, ch, s));
        subj_ha = squeeze(grand_ersp_HA(freq_idx, time_idx, ch, s));

        extracted_data(s, col_idx)   = mean(subj_rp(:));
        extracted_data(s, col_idx+1) = mean(subj_rn(:));
        extracted_data(s, col_idx+2) = mean(subj_la(:));
        extracted_data(s, col_idx+3) = mean(subj_ha(:));

        col_idx = col_idx + 4;
    end
end

T = array2table(extracted_data, 'VariableNames', var_names);
writetable(T, save_excel_path);
fprintf('\n--->success!\n');
