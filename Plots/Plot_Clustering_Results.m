
%% This script loads source data and recreates clustering result figures

clear
close all


% Load source data
fprintf('Loading source data...\n');
load('clustering_source_data.mat');

elecCols = SourceData.elecCols;
Ncls_plot = SourceData.Ncls_plot;
Cn = length(Ncls_plot);

%% Create one large figure with all subjects 


figure('Name', 'All Subjects - Clusters and Burst Density');


all_positions = [];
all_network_labels = [];  % 1=Cue, 2=Target, 3=Both, 0=Neither

for iSubj = 2:12
    subj_data = SourceData.Subjects(iSubj);
    datasetName = subj_data.name;
    subj_label = sprintf('S%d', iSubj);
    col_idx = iSubj - 1;  
  
    optimal_K_idx = subj_data.Optimal_K;
    positions = subj_data.electrode_positions;
    cluster_assignments = subj_data.optimal_clustering.cluster_assignments;
    channel_order = subj_data.optimal_clustering.channel_order;
    n_clusters = Ncls_plot(optimal_K_idx);
    time_vec = subj_data.time_vector;
    
    electrode_network = zeros(size(positions, 1), 1);  % 0=neither, 1=cue, 2=target, 3=both
    
    for iElec = 1:size(positions, 1)
        elec_cluster = cluster_assignments(channel_order == iElec);
        if isempty(elec_cluster)
            continue;
        end
        
        cue_sig = subj_data.burst_density.Cue.cluster(elec_cluster).significant;
        cue_mean = subj_data.burst_density.Cue.cluster(elec_cluster).mean;

        cue_time_window = (time_vec >= 50 & time_vec <= 400);  
        cue_positive = mean(cue_mean(cue_time_window)) > 0;
        cue_network = cue_sig && cue_positive;
        
        target_sig = subj_data.burst_density.Target.cluster(elec_cluster).significant;
        target_mean = subj_data.burst_density.Target.cluster(elec_cluster).mean;

        target_time_window = (time_vec >= 50 & time_vec <= 1000); 
        target_positive = mean(target_mean(target_time_window)) > 0;
        target_network = target_sig && target_positive;
        
        if cue_network && target_network
            electrode_network(iElec) = 3;  % Both
        elseif cue_network
            electrode_network(iElec) = 1;  % Cue only
        elseif target_network
            electrode_network(iElec) = 2;  % Target only
        else
            electrode_network(iElec) = 0;  % Neither
        end
    end
    
    all_positions = [all_positions; positions];
    all_network_labels = [all_network_labels; electrode_network];
    
    %% Row 1: 3D Electrode Scatter Plot
    subplot(4, 11, col_idx)
    hold on
    
    for iClust = 1:n_clusters
        cluster_idx = channel_order(cluster_assignments == iClust);
        scatter3(positions(cluster_idx, 1), ...
                positions(cluster_idx, 2), ...
                positions(cluster_idx, 3), ...
                60, 'filled', ...
                'MarkerFaceColor', elecCols(iClust, :), ...
                'MarkerFaceAlpha', 0.8, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.3);
    end
    
    xlabel('X', 'FontSize', 8)
    ylabel('Y', 'FontSize', 8)
    zlabel('Z', 'FontSize', 8)
    title(sprintf('%s (K=%d)', subj_label, n_clusters), 'FontSize', 9, 'FontWeight', 'bold')
    grid on
    axis equal
    view(45, 30)
    set(gca, 'TickDir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    
    %% Row 2: Cue Epoch Burst Density
    subplot(4, 11, 11 + col_idx)
    hold on
    
    epoch_data = subj_data.burst_density.Cue;
    
    for iClust = 1:n_clusters
        cluster_data = epoch_data.cluster(iClust);
        
        x_fill = [time_vec, fliplr(time_vec)];
        y_fill = [cluster_data.mean + cluster_data.se, ...
                 fliplr(cluster_data.mean - cluster_data.se)];
        
        fill(x_fill, y_fill, elecCols(iClust, :), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        if cluster_data.significant
            plot(time_vec, cluster_data.mean, '-', ...
                'Color', elecCols(iClust, :), 'LineWidth', 1.5);
        else
            plot(time_vec, cluster_data.mean, '-', ...
                'Color', elecCols(iClust, :), 'LineWidth', 0.8);
        end
    end
    
    plot([0 0], ylim, 'k--', 'LineWidth', 0.5)
    plot(xlim, [0 0], 'k--', 'LineWidth', 0.5)
    
    box off
    xlim([-1000 1500])
    if col_idx == 1
        ylabel('Cue (Z-score)', 'FontSize', 9, 'FontWeight', 'bold')
    end
    set(gca, 'TickDir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    set(gca, 'XTickLabel', [])
    
    %% Row 3: Target Epoch Burst Density
    subplot(4, 11, 22 + col_idx)
    hold on
    
    epoch_data = subj_data.burst_density.Target;
    
    for iClust = 1:n_clusters
        cluster_data = epoch_data.cluster(iClust);
        
        x_fill = [time_vec, fliplr(time_vec)];
        y_fill = [cluster_data.mean + cluster_data.se, ...
                 fliplr(cluster_data.mean - cluster_data.se)];
        
        fill(x_fill, y_fill, elecCols(iClust, :), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        if cluster_data.significant
            plot(time_vec, cluster_data.mean, '-', ...
                'Color', elecCols(iClust, :), 'LineWidth', 1.5);
        else
            plot(time_vec, cluster_data.mean, '-', ...
                'Color', elecCols(iClust, :), 'LineWidth', 0.8);
        end
    end
    
    plot([0 0], ylim, 'k--', 'LineWidth', 0.5)
    plot(xlim, [0 0], 'k--', 'LineWidth', 0.5)
    
    
    box off
    xlim([-1000 1500])
    if col_idx == 1
        ylabel('Target (Z-score)', 'FontSize', 9, 'FontWeight', 'bold')
    end
    xlabel('Time (ms)', 'FontSize', 8)
    set(gca, 'TickDir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
end
set(gcf,'Position',[80   27   1570   1250])
%% Row 4: Combined Network Plot (spanning all columns)
figure
hold on
mesh = ft_read_headshape('surface_white_both.mat'); 
ft_plot_mesh(mesh, 'edgecolor','none', 'facealpha',0.4,'facecolor',[.9 .9 .9]);
view(0, 90); 
camlight; material dull; lighting gouraud;


purple = [0.6, 0.2, 0.8];  % Cue network
green = [0.2, 0.7, 0.2];   % Target network
yellow = [0.9, 0.7, 0.2];  % Both networks
gray = [0.4, 0.4, 0.4];    % Neither

for net_type = 0:3
    idx = find(all_network_labels == net_type);
    if isempty(idx)
        continue;
    end
    
    switch net_type
        case 0 % Neither
            color = gray;
            marker_size = 70;
            alpha = 0.9;
            label_str = 'Neither';
        case 1 % Cue only
            color = purple;
            marker_size = 70;
            alpha = 0.9;
            label_str = 'Cue Network';
        case 2 % Target only
            color = green;
            marker_size = 70;
            alpha = 0.9;
            label_str = 'Target Network';
        case 3 % Both
            color = yellow;
            marker_size = 60;
            alpha = 0.9;
            label_str = 'Both Networks';
    end
    
    scatter3(all_positions(idx, 1), ...
            all_positions(idx, 2), ...
            all_positions(idx, 3), ...
            marker_size, 'filled', ...
            'MarkerFaceColor', color, ...
            'MarkerFaceAlpha', alpha, ...
            'MarkerEdgeColor', 'none', ...
            'LineWidth', 0.5, ...
            'DisplayName', label_str);
end

xlabel('X (mm)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Y (mm)', 'FontSize', 11, 'FontWeight', 'bold')
zlabel('Z (mm)', 'FontSize', 11, 'FontWeight', 'bold')
title('Cue and Target Networks Across All Subjects', 'FontSize', 12, 'FontWeight', 'bold')
legend('Location', 'eastoutside', 'Box', 'off', 'FontSize', 10)
grid on
axis equal
set(gca, 'TickDir', 'out', 'FontName', 'Helvetica', 'FontSize', 10)

