% This script loads the source data and creates:
%   1. Lead-lag analysis plot
%   2. Burst-triggered HFB time-series plots (shaded error bars)


clear
close all


% Load source data
load( 'leadlag_source_data.mat');
fprintf('Data loaded successfully.\n\n');



figure('Name', 'Lead-Lag Analysis for Both Tasks', 'Position', [100 400 1000 500]);

% Task 1 Plot (Subjects 1-7)
subplot(1,2,1)
title('Task 1 (Subjects 1-7)', 'FontSize', 12, 'FontWeight', 'bold')
hold on
grid on

task1_subjects = SourceData.info.task1_subjects;

for idx = 1:length(task1_subjects)
    iSubj = task1_subjects(idx);
    subj = SourceData.subjects{iSubj};
    
    Mn1 = subj.stats.Att_to_Tar_median;
    SE1 = subj.stats.Att_to_Tar_SE;
    Mn2 = subj.stats.Tar_to_Att_median;
    SE2 = subj.stats.Tar_to_Att_SE;
    

    plot([Mn1-SE1, Mn1+SE1], [idx, idx], 'LineWidth', 2, 'Color', 'r')
    plot([Mn2-SE2, Mn2+SE2], [idx, idx], 'LineWidth', 2, 'Color', 'b')
    
    plot(Mn1, idx, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 6)
    plot(Mn2, idx, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 6)
end

xlabel('Peak Time (ms)', 'FontSize', 11)
ylabel('Subject', 'FontSize', 11)
legend({'Att→Tar', 'Tar→Att'}, 'Location', 'best', 'FontSize', 10)
xlim([-250 250])
ylim([0.5 7.5])
set(gca, 'YDir', 'reverse')
set(gca, 'YTick', 1:7)

plot([0 0], ylim, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off')

% Task 2 Plot (Subjects 8-12)
subplot(1,2,2)
title('Task 2 (Subjects 8-12)', 'FontSize', 12, 'FontWeight', 'bold')
hold on
grid on

task2_subjects = SourceData.info.task2_subjects;

for idx = 1:length(task2_subjects)
    iSubj = task2_subjects(idx);
    subj = SourceData.subjects{iSubj};
    
    Mn1 = subj.stats.Att_to_Tar_median;
    SE1 = subj.stats.Att_to_Tar_SE;
    Mn2 = subj.stats.Tar_to_Att_median;
    SE2 = subj.stats.Tar_to_Att_SE;
    
    plot([Mn1-SE1, Mn1+SE1], [idx, idx], 'LineWidth', 2, 'Color', 'r')
    plot([Mn2-SE2, Mn2+SE2], [idx, idx], 'LineWidth', 2, 'Color', 'b')
    
    plot(Mn1, idx, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 6)
    plot(Mn2, idx, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 6)
end

xlabel('Peak Time (ms)', 'FontSize', 11)
ylabel('Subject', 'FontSize', 11)
legend({'Att→Tar', 'Tar→Att'}, 'Location', 'best', 'FontSize', 10)
xlim([-250 250])
ylim([0.5 5.5])
set(gca, 'YDir', 'reverse')
set(gca, 'YTick', 1:5)

plot([0 0], ylim, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off')


% FIGURE 2: BURST-TRIGGERED HFB - TASK 1 (Shaded Error Bars)

figure('Name', 'Task 1 - Burst-Triggered HFB (Subjects 1-7)', 'Position', [100 100 1400 800]);

for idx = 1:7
    iSubj = idx;
    subj = SourceData.subjects{iSubj};
    
    subplot(2, 4, idx)
    hold on
    grid on
    
    t_vec = subj.timeseries.time_vector;
    
    if ~isempty(subj.timeseries.Att_to_Tar_mean)
        att_tar_mean = subj.timeseries.Att_to_Tar_mean; 
        att_tar_sem = subj.timeseries.Att_to_Tar_sem;   

        if size(att_tar_mean, 1) > 1
            att_tar_mean = att_tar_mean';
            att_tar_sem = att_tar_sem';
        end
        
        fill([t_vec, fliplr(t_vec)], ...
             [att_tar_mean + att_tar_sem, fliplr(att_tar_mean - att_tar_sem)], ...
             [0.8, 0.2, 0.2], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        plot(t_vec, att_tar_mean, 'Color', [0.8, 0.2, 0.2], 'LineWidth', 2);
    end
    
    if ~isempty(subj.timeseries.Tar_to_Att_mean)
        tar_att_mean = subj.timeseries.Tar_to_Att_mean; 
        tar_att_sem = subj.timeseries.Tar_to_Att_sem;   

        if size(tar_att_mean, 1) > 1
            tar_att_mean = tar_att_mean';
            tar_att_sem = tar_att_sem';
        end
        
        fill([t_vec, fliplr(t_vec)], ...
             [tar_att_mean + tar_att_sem, fliplr(tar_att_mean - tar_att_sem)], ...
             [0.2, 0.2, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        plot(t_vec, tar_att_mean, 'Color', [0.2, 0.2, 0.8], 'LineWidth', 2);
    end
    
    xlabel('Time to burst (ms)', 'FontSize', 10)
    ylabel('Normalized HFB', 'FontSize', 10)
    title(sprintf('S%d: %s', iSubj, subj.name), 'FontSize', 11, 'FontWeight', 'bold')
    xlim([-500 500])
    
    plot([0 0], ylim, 'k--', 'LineWidth', 0.5)
    plot(xlim, [0 0], 'k--', 'LineWidth', 0.5)
    
    if idx == 1
        legend({'Att→Tar', '', 'Tar→Att', ''}, 'Location', 'best', 'FontSize', 9)
    end
end


% FIGURE 3: BURST-TRIGGERED HFB - TASK 2 (Shaded Error Bars)

figure('Name', 'Task 2 - Burst-Triggered HFB (Subjects 8-12)', 'Position', [100 100 1400 500]);

for idx = 1:5
    iSubj = idx + 7;
    subj = SourceData.subjects{iSubj};
    
    subplot(2, 3, idx)
    hold on
    grid on
    
    t_vec = subj.timeseries.time_vector; 
    
    if ~isempty(subj.timeseries.Att_to_Tar_mean)
        att_tar_mean = subj.timeseries.Att_to_Tar_mean; 
        att_tar_sem = subj.timeseries.Att_to_Tar_sem;   
        
        if size(att_tar_mean, 1) > 1
            att_tar_mean = att_tar_mean';
            att_tar_sem = att_tar_sem';
        end
        
        fill([t_vec, fliplr(t_vec)], ...
             [att_tar_mean + att_tar_sem, fliplr(att_tar_mean - att_tar_sem)], ...
             [0.8, 0.2, 0.2], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        plot(t_vec, att_tar_mean, 'Color', [0.8, 0.2, 0.2], 'LineWidth', 2);
    end
    
    if ~isempty(subj.timeseries.Tar_to_Att_mean)
        tar_att_mean = subj.timeseries.Tar_to_Att_mean; 
        tar_att_sem = subj.timeseries.Tar_to_Att_sem;   
        

        if size(tar_att_mean, 1) > 1
            tar_att_mean = tar_att_mean';
            tar_att_sem = tar_att_sem';
        end
        

        fill([t_vec, fliplr(t_vec)], ...
             [tar_att_mean + tar_att_sem, fliplr(tar_att_mean - tar_att_sem)], ...
             [0.2, 0.2, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        

        plot(t_vec, tar_att_mean, 'Color', [0.2, 0.2, 0.8], 'LineWidth', 2);
    end
    

    xlabel('Time to burst (ms)', 'FontSize', 10)
    ylabel('Normalized HFB', 'FontSize', 10)
    title(sprintf('S%d: %s', iSubj, subj.name), 'FontSize', 11, 'FontWeight', 'bold')
    xlim([-500 500])
    

    plot([0 0], ylim, 'k--', 'LineWidth', 0.5)
    plot(xlim, [0 0], 'k--', 'LineWidth', 0.5)
    

    if idx == 1
        legend({'Att→Tar', '', 'Tar→Att', ''}, 'Location', 'best', 'FontSize', 9)
    end
end

fprintf('\n=== Summary Statistics ===\n\n');

fprintf('Task 1 (Subjects 1-7):\n');
for idx = 1:length(task1_subjects)
    iSubj = task1_subjects(idx);
    subj = SourceData.subjects{iSubj};
    fprintf('  Subject %d (%s):\n', iSubj, subj.name);
    fprintf('    Att channels: %d, Tar channels: %d\n', ...
            sum(subj.Att_channels), sum(subj.Tar_channels));
    fprintf('    Att→Tar: %.1f ± %.1f ms\n', subj.stats.Att_to_Tar_median, subj.stats.Att_to_Tar_SE);
    fprintf('    Tar→Att: %.1f ± %.1f ms\n', subj.stats.Tar_to_Att_median, subj.stats.Tar_to_Att_SE);
    if ~isnan(subj.stats.ranksum_p)
        fprintf('    p-value: %.4f\n', subj.stats.ranksum_p);
    end
end

fprintf('\nTask 2 (Subjects 8-12):\n');
for idx = 1:length(task2_subjects)
    iSubj = task2_subjects(idx);
    subj = SourceData.subjects{iSubj};
    fprintf('  Subject %d (%s):\n', iSubj, subj.name);
    fprintf('    Att channels: %d, Tar channels: %d\n', ...
            sum(subj.Att_channels), sum(subj.Tar_channels));
    fprintf('    Att→Tar: %.1f ± %.1f ms\n', subj.stats.Att_to_Tar_median, subj.stats.Att_to_Tar_SE);
    fprintf('    Tar→Att: %.1f ± %.1f ms\n', subj.stats.Tar_to_Att_median, subj.stats.Tar_to_Att_SE);
    if ~isnan(subj.stats.ranksum_p)
        fprintf('    p-value: %.4f\n', subj.stats.ranksum_p);
    end
end

