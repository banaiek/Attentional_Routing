%% Plot Mean btHFA, PC1, and PC2 Time Courses

clear
close all


% Load source data
fprintf('Loading PC time course source data...\n');
load([ 'SourceData_PC_TimeCourses.mat']);


time = SourceData_TimeCourses.time;


%% Plot 2: All Subjects Individual Plots
figure('Position', [100, 100, 1600, 1000], 'Color', 'w');

num_subjects = length(SourceData_TimeCourses.subjects);
rows = 4;
cols = 3;

for iSubj = 1:num_subjects
    
    subplot(rows, cols, iSubj);
    
    btHFA_norm = SourceData_TimeCourses.mean_btHFA{iSubj} / ...
                 max(abs(SourceData_TimeCourses.mean_btHFA{iSubj}));
    PC1_norm = SourceData_TimeCourses.PC1{iSubj} / ...
               max(abs(SourceData_TimeCourses.PC1{iSubj}));
    PC2_norm = SourceData_TimeCourses.PC2{iSubj} / ...
               max(abs(SourceData_TimeCourses.PC2{iSubj}));
    
    plot(time, btHFA_norm, 'k-', 'LineWidth', 1.5); hold on;
    plot(time, PC1_norm, 'r-', 'LineWidth', 1.5);
    plot(time, PC2_norm, 'b-', 'LineWidth', 1.5);
    plot([0 0], ylim, 'k--', 'LineWidth', 0.5);
    hold off;
    
    title(sprintf('%s (PC1:%.1f%%, PC2:%.1f%%)', ...
                  SourceData_TimeCourses.subjects{iSubj}, ...
                  SourceData_TimeCourses.PC1_variance_explained(iSubj), ...
                  SourceData_TimeCourses.PC2_variance_explained(iSubj)), ...
          'FontSize', 10);
    
    if iSubj > (rows-1)*cols
        xlabel('Time (ms)', 'FontSize', 9);
    end
    if mod(iSubj-1, cols) == 0
        ylabel('Normalized Amplitude', 'FontSize', 9);
    end
    
    grid on;
    xlim([time(1) time(end)]);
    set(gca, 'FontSize', 8);
    
    if iSubj == 1
        legend({'Mean btHFA', 'PC1', 'PC2'}, 'Location', 'best', 'FontSize', 8);
    end
end





