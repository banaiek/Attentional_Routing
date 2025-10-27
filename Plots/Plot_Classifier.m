clear
close all

%% Load source data
load( 'Classifier_SourceData.mat');

time_axis = SourceData.metadata.classifier_Timepoints + SourceData.metadata.classifier_Win/2;

response_types = {'CueResponsive', 'CueNonResponsive'};
response_labels = {'Cue Responsive', 'Cue Non-Responsive'};


figure('Position', [100 100 400 500]);

for iResp = 1:2
    resp_name = response_types{iResp};
    
    accuracy = SourceData.(resp_name).Ts_Accuracy;
    se = SourceData.(resp_name).seAccuracy;
    pValues = SourceData.(resp_name).adjPValues;
    
    if iResp == 1
        color = [0.4 0.4470 0.910]; 
    else
        color = [0.500 0.4250 0.4980]; 
    end
    
    h = shadedErrorBar(time_axis, accuracy, se, ...
                  {'Color', color, 'LineWidth', 2.5, ...
                               'DisplayName', response_labels{iResp}});
    

    set(h.edge, 'HandleVisibility', 'off');
    hold on;
    
    hVal = double(pValues < 0.05);
     hVal(time_axis <= 0) = 0; 
     hVal(hVal == 0) = nan;
    plot(time_axis, accuracy .* hVal, 'LineWidth', 3, ...
         'Color', 'r', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, ...
         'HandleVisibility', 'off');
end


chanceLevel = SourceData.CueResponsive.chanceLevel;
plot(xlim, [chanceLevel chanceLevel], 'LineWidth', 1.5, 'Color', 'k', ...
     'LineStyle', '--', 'DisplayName', 'Chance Level');

plot([0 0], ylim, 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', ...
     'DisplayName', 'Cue Onset');

xlim([-500 750]);
xlabel('Time (ms)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Classification Accuracy', 'FontSize', 13, 'FontWeight', 'bold');
title('Task 1: Cue Responsive vs Non-Responsive Channels', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1.5);
box off;
grid on;
legend('Location', 'best', 'FontSize', 11);

