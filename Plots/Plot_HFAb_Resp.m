%% Generate HFAb resp Figures from Source Data
clear
close all

%% Load source data

fprintf('Loading source data...\n');
load( 'HFA_Resp_source_data.mat');

time_vector = source_data.time_vector;
epoch_labels = source_data.epoch_labels;


Col_Laterality = [.2 .7 .2; .7 .2 .2];  
Col_Validity = [.2 .2 .7; .7 .2 .2];    

XLbl = {'Time to Cue (ms)', 'Time to Target (ms)', 'Time to Response (ms)'};
yLim_Egly = [-.6 3.5];
yLim_STN = [-.6 2];

%% Figure 1: Laterality (CL vs IL)
fprintf('Generating laterality plots...\n');

figure('Position', [168 338 1411 540]);


for iEpoch = 1:length(epoch_labels)
    subplot(2, 3, iEpoch)
    hold on
    

    for iCond = 1:length(source_data.laterality.labels)
        Mn = source_data.laterality.egly.mean{iCond, iEpoch};
        SE = source_data.laterality.egly.sem{iCond, iEpoch};
        

        s(iCond) = shadedErrorBar(time_vector, Mn, SE, ...
            {'color', Col_Laterality(iCond,:), 'LineWidth', 2}, 1);
    end
    

    box off
    set(gca, 'TickDir', 'out', 'XLim', [-1000 1500], 'YLim', yLim_Egly)
    
    plot([0 0], get(gca, 'ylim'), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    plot(get(gca, 'xlim'), [0 0], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    
    xlabel(XLbl{iEpoch})
    ylabel('Burst Density (Z-score)')
    title(sprintf('Laterality\nEgly'))
    legend([s(1).patch, s(2).patch], source_data.laterality.labels, 'Box', 'off')
end

for iEpoch = 1:length(epoch_labels)
    subplot(2, 3, iEpoch + 3)
    hold on
    
    for iCond = 1:length(source_data.laterality.labels)
        Mn = source_data.laterality.stn.mean{iCond, iEpoch};
        SE = source_data.laterality.stn.sem{iCond, iEpoch};
        
        s(iCond) = shadedErrorBar(time_vector, Mn, SE, ...
            {'color', Col_Laterality(iCond,:), 'LineWidth', 2}, 1);
    end
    
    box off
    set(gca, 'TickDir', 'out', 'XLim', [-1000 1500], 'YLim', yLim_STN)
    
    plot([0 0], get(gca, 'ylim'), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    plot(get(gca, 'xlim'), [0 0], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    
    xlabel(XLbl{iEpoch})
    ylabel('Burst Density (Z-score)')
    title(sprintf('Laterality\nSTN'))
    legend([s(1).patch, s(2).patch], source_data.laterality.labels, 'Box', 'off')
end



%% Figure 2: Validity (Valid vs Invalid)
fprintf('Generating validity plots...\n');

figure('Position', [168 338 1411 540]);

for iEpoch = 1:length(epoch_labels)
    subplot(2, 3, iEpoch)
    hold on
    
    for iCond = 1:length(source_data.validity.labels)
        Mn = source_data.validity.egly.mean{iCond, iEpoch};
        SE = source_data.validity.egly.sem{iCond, iEpoch};
        
        s(iCond) = shadedErrorBar(time_vector, Mn, SE, ...
            {'color', Col_Validity(iCond,:), 'LineWidth', 2}, 1);
    end
    
    box off
    set(gca, 'TickDir', 'out', 'XLim', [-1000 1500], 'YLim', yLim_Egly)
    
    plot([0 0], get(gca, 'ylim'), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    plot(get(gca, 'xlim'), [0 0], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    
    xlabel(XLbl{iEpoch})
    ylabel('Burst Density (Z-score)')
    title(sprintf('Validity\nEgly'))
    legend([s(1).patch, s(2).patch], source_data.validity.labels, 'Box', 'off')
end

for iEpoch = 1:length(epoch_labels)
    subplot(2, 3, iEpoch + 3)
    hold on
    
    for iCond = 1:length(source_data.validity.labels)
        Mn = source_data.validity.stn.mean{iCond, iEpoch};
        SE = source_data.validity.stn.sem{iCond, iEpoch};
        
        s(iCond) = shadedErrorBar(time_vector, Mn, SE, ...
            {'color', Col_Validity(iCond,:), 'LineWidth', 2}, 1);
    end
    
    box off
    set(gca, 'TickDir', 'out', 'XLim', [-1000 1500], 'YLim', yLim_STN)
    
    plot([0 0], get(gca, 'ylim'), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    plot(get(gca, 'xlim'), [0 0], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    
    xlabel(XLbl{iEpoch})
    ylabel('Burst Density (Z-score)')
    title(sprintf('Validity\nSTN'))
    legend([s(1).patch, s(2).patch], source_data.validity.labels, 'Box', 'off')
end

