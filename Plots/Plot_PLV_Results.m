%% Plot PLV Results from Source Data


clear
clc



% Load source data
load('PLV_source_data.mat');

Steps = source_data.metadata.time_steps;
freqs = source_data.metadata.frequencies;
TOI = source_data.metadata.TOI;
FOI = source_data.metadata.FOI;
nSubjects = source_data.metadata.n_subjects;
task_assignment = source_data.metadata.task_assignment;

PLV_Cue = source_data.PLV_Cue;
PLV_Target = source_data.PLV_Target;
Rand_Max_Cue = source_data.Rand_Max_Cue;
Rand_Max_Target = source_data.Rand_Max_Target;
Rand_Min_Cue = source_data.Rand_Min_Cue;
Rand_Min_Target = source_data.Rand_Min_Target;

theta_alpha_band = find(freqs >= 4 & freqs <= 12);
beta_band = find(freqs >= 15 & freqs <= 25);

XLbl = {'Time to Cue (ms)', 'Time to Target (ms)'};
task_names = {'Task 1', 'Task 2'};

Sig_DeCoupled_Cue = cellfun(@(x,y) x < prctile(y, 1, [1]), ...
    PLV_Cue, Rand_Min_Cue, 'UniformOutput', false);
Sig_DeCoupled_Target = cellfun(@(x,y) x < prctile(y, 1, [1]), ...
    PLV_Target, Rand_Min_Target, 'UniformOutput', false);



%% Figure 1: Average PLV by Task
fig2 = figure('Position', [100, 100, 1200, 800]);
fig2.Name = 'PLV Time-Frequency - Task Average';

% Task 1 subjects (1-7)
task1_subjects = find(task_assignment == 1);
% Task 2 subjects (8-12)
task2_subjects = find(task_assignment == 2);

for iEpoch = 1:2  % Cue and Target
    if iEpoch == 1
        PLV_data = PLV_Cue;
    else
        PLV_data = PLV_Target;
    end
    
    % Task 1
    subplot(2, 2, iEpoch);
 
    PLV_task1_cell = cellfun(@(x) nanmean(x, 3), PLV_data(task1_subjects), 'UniformOutput', false);
    PLV_task1 = cat(3, PLV_task1_cell{:});
    PLV_task1_avg = nanmean(PLV_task1, 3);  
    PLV_task1_norm = (PLV_task1_avg - nanmean(PLV_task1_avg(TOI, FOI), [1,2])) ./ ...
        std(PLV_task1_avg(TOI, FOI), [], [1,2], 'omitmissing');
    
    imagesc(Steps(TOI), freqs(FOI), PLV_task1_norm(TOI, FOI)');
    axis xy;
    colormap jet;
    hold on;
    plot([0 0], get(gca, 'ylim'), 'k--', 'LineWidth', 1.5);
    xlabel(XLbl{iEpoch});
    ylabel('Frequency (Hz)');
    title(sprintf('Task 1 (n=%d) - %s', length(task1_subjects), XLbl{iEpoch}));
    set(gca, 'TickDir', 'out', 'XLim', [-1250 1250], 'CLim', [-1 2.5]);
    colorbar;
    
    % Task 2
    subplot(2, 2, iEpoch + 2);

    PLV_task2_cell = cellfun(@(x) nanmean(x, 3), PLV_data(task2_subjects), 'UniformOutput', false);
    PLV_task2 = cat(3, PLV_task2_cell{:});
    PLV_task2_avg = nanmean(PLV_task2, 3);  
    PLV_task2_norm = (PLV_task2_avg - nanmean(PLV_task2_avg(TOI, FOI), [1,2])) ./ ...
        std(PLV_task2_avg(TOI, FOI), [], [1,2], 'omitmissing');
    
    imagesc(Steps(TOI), freqs(FOI), PLV_task2_norm(TOI, FOI)');
    axis xy;
    colormap jet;
    hold on;
    plot([0 0], get(gca, 'ylim'), 'k--', 'LineWidth', 1.5);
    xlabel(XLbl{iEpoch});
    ylabel('Frequency (Hz)');
    title(sprintf('Task 2 (n=%d) - %s', length(task2_subjects), XLbl{iEpoch}));
    set(gca, 'TickDir', 'out', 'XLim', [-1250 1250], 'CLim', [-1 2.5]);
    colorbar;
end



%% Figure 2: Significant Decoupling - Theta/Alpha and Beta Bands (Line Plots Only)
fig3 = figure('Position', [100, 100, 1200, 600]);
fig3.Name = 'Significant PLV Decoupling - Band Averages';

for iEpoch = 1:2
    if iEpoch == 1
        Sig_data = Sig_DeCoupled_Cue;
        epoch_name = 'Cue';
    else
        Sig_data = Sig_DeCoupled_Target;
        epoch_name = 'Target';
    end

    subplot(2, 2, (iEpoch-1)*2 + 1);

    Sig_task1_cell = cellfun(@(x) nanmean(x, 3), Sig_data(task1_subjects), 'UniformOutput', false);
    Sig_task1 = cat(3, Sig_task1_cell{:});
    Sig_task1_avg = squeeze(nanmean(Sig_task1(TOI, :, :), 3));
    

    baseline_idx = find(Steps(TOI) < 0);
    post_zero_idx = find(Steps(TOI) >= 0);
    baseline_t1 = nanmean(Sig_task1_avg(baseline_idx, :), 1);
    Sig_task1_corrected = Sig_task1_avg - baseline_t1 + 0.05;
    

    Sig_theta_alpha = squeeze(nanmean(Sig_task1_corrected(:, theta_alpha_band), 2));
    Sig_beta = squeeze(nanmean(Sig_task1_corrected(:, beta_band), 2));
    

    Sig_theta_alpha_smooth = smooth(Sig_theta_alpha, 5);
    Sig_beta_smooth = smooth(Sig_beta, 5);
    

    baseline_theta_alpha = nanmean(Sig_theta_alpha_smooth(baseline_idx));
    baseline_beta = nanmean(Sig_beta_smooth(baseline_idx));
    std_theta_alpha = nanstd(Sig_theta_alpha_smooth(baseline_idx));
    std_beta = nanstd(Sig_beta_smooth(baseline_idx));
    thresh_theta_alpha = baseline_theta_alpha + 2*std_theta_alpha;
    thresh_beta = baseline_beta + 2*std_beta;
    
    hold on;
    
    for i = 1:length(TOI)
        if i > 1 && Steps(TOI(i)) >= 0 && Sig_theta_alpha(i) > thresh_theta_alpha
            plot(Steps(TOI(i-1:i)), Sig_theta_alpha_smooth(i-1:i), 'color',[.2 .8 .3], 'LineWidth', 4);
        else

            if i > 1
                plot(Steps(TOI(i-1:i)), Sig_theta_alpha_smooth(i-1:i), 'color',[.2 .8 .3], 'LineWidth', 2);
            end
        end
    end
    

    for i = 1:length(TOI)
        if i > 1 && Steps(TOI(i)) >= 0 && Sig_beta(i) > thresh_beta
            plot(Steps(TOI(i-1:i)), Sig_beta_smooth(i-1:i), 'color',[.2 .3 .8], 'LineWidth', 4);
        else

            if i > 1
                plot(Steps(TOI(i-1:i)), Sig_beta_smooth(i-1:i), 'color',[.2 .3 .8], 'LineWidth', 2);
            end
        end
    end
    
    plot([0 0], get(gca, 'ylim'), 'k--', 'LineWidth', 1.5);
    yline(0.05, 'k:', 'LineWidth', 1.5, 'Label', 'Baseline (0.05)');
    xlabel(sprintf('Time to %s (ms)', epoch_name));
    ylabel('Proportion Significant');
    title(sprintf('Task 1 - %s - Band Average', epoch_name));
    legend('Theta/Alpha (4-12 Hz)', 'Beta (15-25 Hz)', 'Location', 'best');
    set(gca, 'TickDir', 'out', 'XLim', [-1250 1250]);
    grid on;
    

    subplot(2, 2, (iEpoch-1)*2 + 2);

    Sig_task2_cell = cellfun(@(x) nanmean(x, 3), Sig_data(task2_subjects), 'UniformOutput', false);
    Sig_task2 = cat(3, Sig_task2_cell{:});
    Sig_task2_avg = squeeze(nanmean(Sig_task2(TOI, :, :), 3));
    

    baseline_t2 = nanmean(Sig_task2_avg(baseline_idx, :), 1);
    Sig_task2_corrected = Sig_task2_avg - baseline_t2 + 0.05;
    

    Sig_theta_alpha_t2 = squeeze(nanmean(Sig_task2_corrected(:, theta_alpha_band), 2));
    Sig_beta_t2 = squeeze(nanmean(Sig_task2_corrected(:, beta_band), 2));
    

    Sig_theta_alpha_t2_smooth = smooth(Sig_theta_alpha_t2, 5);
    Sig_beta_t2_smooth = smooth(Sig_beta_t2, 5);
    

    baseline_theta_alpha_t2 = nanmean(Sig_theta_alpha_t2_smooth(baseline_idx));
    baseline_beta_t2 = nanmean(Sig_beta_t2_smooth(baseline_idx));
    std_theta_alpha_t2 = nanstd(Sig_theta_alpha_t2_smooth(baseline_idx));
    std_beta_t2 = nanstd(Sig_beta_t2_smooth(baseline_idx));
    thresh_theta_alpha_t2 = baseline_theta_alpha_t2 + 2*std_theta_alpha_t2;
    thresh_beta_t2 = baseline_beta_t2 + 2*std_beta_t2;
    
    hold on;
    
    for i = 1:length(TOI)
        if i > 1 && Steps(TOI(i)) >= 0 && Sig_theta_alpha_t2(i) > thresh_theta_alpha_t2
            plot(Steps(TOI(i-1:i)), Sig_theta_alpha_t2_smooth(i-1:i), 'color',[.2 .8 .3], 'LineWidth', 4);
        else

            if i > 1
                plot(Steps(TOI(i-1:i)), Sig_theta_alpha_t2_smooth(i-1:i), 'color',[.2 .8 .3], 'LineWidth', 2);
            end
        end
    end
    
    for i = 1:length(TOI)
        if i > 1 && Steps(TOI(i)) >= 0 && Sig_beta_t2(i) > thresh_beta_t2
            plot(Steps(TOI(i-1:i)), Sig_beta_t2_smooth(i-1:i), 'color',[.2 .3 .8], 'LineWidth', 4);
        else

            if i > 1
                plot(Steps(TOI(i-1:i)), Sig_beta_t2_smooth(i-1:i), 'color',[.2 .3 .8], 'LineWidth', 2);
            end
        end
    end
    
    plot([0 0], get(gca, 'ylim'), 'k--', 'LineWidth', 1.5);
    yline(0.05, 'k:', 'LineWidth', 1.5, 'Label', 'Baseline (0.05)');
    xlabel(sprintf('Time to %s (ms)', epoch_name));
    ylabel('Proportion Significant');
    title(sprintf('Task 2 - %s - Band Average', epoch_name));
    legend('Theta/Alpha (4-12 Hz)', 'Beta (15-25 Hz)', 'Location', 'best');
    set(gca, 'TickDir', 'out', 'XLim', [-1250 1250]);
    grid on;
end

