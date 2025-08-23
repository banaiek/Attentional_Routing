function RESClust = perform_sampling_clustering(Clust_vec, proportion, num_samples, Ncls, Corr_type)
% two-level sampling-based clustering with consensus
% Inputs:
%   Clust_vec - clustering matrix (channels x channels)
%   proportion - proportion of data to subsample (0-1)
%   num_samples - number of iterations
%   Ncls - vector of cluster numbers to test
%   Corr_type - distance metric for kmeans
% Outputs:
%   RESClust - structure containing all clustering results

    Num_Channels = size(Clust_vec, 1);
    
    if isempty(gcp('nocreate'))
        parpool;
    end
    
    stream = RandStream('mlfg6331_64');
    options = statset('UseParallel', 1, 'UseSubstreams', 1, 'Streams', stream);
    
    % Calculate sample size
    sample_size = max(floor(Num_Channels * proportion), max(Ncls) * 5);
    sample_size = min(sample_size, Num_Channels);
    
    All_ClustRes = cell(length(Ncls), 1);
    P_clust = zeros(length(Ncls), Num_Channels, Num_Channels);
    
    % First level: clustering on original data
    fprintf('  First level clustering...\n');
    for iClust = 1:length(Ncls)
        fprintf('    Testing k=%d clusters\n', Ncls(iClust));
        
        % Generate subsamples
        [subsamples, pair_count] = sampleCount(Num_Channels, sample_size, num_samples);
        
        ClustRes = perform_parallel_clustering(Clust_vec, subsamples, Ncls(iClust), ...
                                              options, Corr_type, num_samples);
        
        All_ClustRes{iClust} = ClustRes;
        P_clust(iClust, :, :) = squeeze(nansum(ClustRes)) ./ pair_count;
    end
    
    % Second level: clustering on consensus matrices
    fprintf('  Second level clustering...\n');
    All_ClustRes_vec = cell(length(Ncls), 1);
    P_clust_vecP = zeros(length(Ncls), Num_Channels, Num_Channels);
    
    for iClust = 1:length(Ncls)
        fprintf('    Testing k=%d clusters\n', Ncls(iClust));
        
        % Generate subsamples
        [subsamples, pair_count] = sampleCount(Num_Channels, sample_size, num_samples);
        
        % Use consensus matrix from first round
        P_clust_Vec = squeeze(P_clust(iClust, :, :));
        P_clust_Vec(isinf(P_clust_Vec)) = 1;
        
        ClustRes_P_vec = perform_parallel_clustering(P_clust_Vec, subsamples, ...
                                                    Ncls(iClust), options, ...
                                                    'sqeuclidean', num_samples);
        
        All_ClustRes_vec{iClust} = ClustRes_P_vec;
        P_clust_vecP(iClust, :, :) = squeeze(nansum(ClustRes_P_vec)) ./ pair_count;
    end
    
    % Order clusters by stability
    [Ordered_Clust, Ordered_Inds, ConfMat, ConfOut] = order_clusters_by_stability(P_clust_vecP, ...
                                                                                  Ncls, options, ...
                                                                                  Corr_type);
    
    RESClust.All_ClustRes = All_ClustRes;
    RESClust.All_ClustRes_vec = All_ClustRes_vec;
    RESClust.P_clust = P_clust;
    RESClust.P_clust_vecP = P_clust_vecP;
    RESClust.Ordered_Clust = Ordered_Clust;
    RESClust.Ordered_Inds = Ordered_Inds;
    RESClust.ConfMat = ConfMat;
    RESClust.Out = ConfOut;
end

function ClustRes = perform_parallel_clustering(data_matrix, subsamples, n_clusters, ...
                                               options, distance_type, num_samples)
% Perform parallel k-means clustering on subsamples
    ClustRes = nan(num_samples, size(data_matrix, 1), size(data_matrix, 2));
    tempRes = cell(num_samples, 1);
    
    parfor iC = 1:num_samples
        indC = sort(subsamples(iC, :));
        P_sub = data_matrix(indC, indC);
        
        % K-means clustering
        [sinds, ~, ~, ~] = kmeans(P_sub, n_clusters, ...
                                  'Options', options, ...
                                  'Replicate', 1000, ...
                                  'MaxIter', 10000, ...
                                  'distance', distance_type);
        
        M_ind = sinds' == sinds;
        tempV = nan(size(data_matrix));
        tempV(indC, indC) = M_ind;
        tempRes{iC} = tempV;
    end
    
    for iC = 1:num_samples
        ClustRes(iC, :, :) = tempRes{iC};
    end
end

function [Ordered_Clust, Ordered_Inds, ConfMat, ConfOut] = order_clusters_by_stability(P_clust_vecP, ...
                                                                                      Ncls, options, ...
                                                                                      Corr_type)
% Order clusters by their stability scores
    Ordered_Clust = zeros(length(Ncls), size(P_clust_vecP, 2));
    Ordered_Inds = zeros(length(Ncls), size(P_clust_vecP, 2));
    ConfMat = cell(length(Ncls), 1);
    ConfOut = cell(length(Ncls), 1);
    
    for iC = 1:length(Ncls)
        mP_clust = squeeze(P_clust_vecP(iC, :, :));
        mP_clust(isinf(mP_clust)) = 1;
        
        % Final clustering on consensus matrix
        [sinds, ~, ~, ~] = kmeans(mP_clust, Ncls(iC), ...
                                  'Options', options, ...
                                  'Replicate', 1000, ...
                                  'MaxIter', 10000, ...
                                  'distance', Corr_type);
        
        % Calculate confusion matrix
        TempConf = zeros(Ncls(iC));
        for i = 1:Ncls(iC)
            for j = 1:Ncls(iC)
                TempConf(i, j) = nanmean(mP_clust(sinds == i, sinds == j), 'all');
            end
        end
        
        % Calculate cluster stability
        Clust_Stability = diag(TempConf) ./ nanmean(TempConf, 2);
        [~, sorted_cls] = sort(Clust_Stability, 'descend');
        
        % Reorder clusters by stability
        temp_sinds = sinds;
        for i = 1:Ncls(iC)
            sinds(temp_sinds == sorted_cls(i)) = i;
        end
        
        % Create ordered indices
        sOrd = [];
        s_ids = [];
        for i = 1:Ncls(iC)
            tOrd = find(sinds == i);
            sOrd = [sOrd; tOrd];
            s_ids = [s_ids; i * ones(size(tOrd))];
        end
        
        Ordered_Clust(iC, :) = sOrd';
        Ordered_Inds(iC, :) = s_ids';
        
        % Final confusion matrix with reordered clusters
        for i = 1:Ncls(iC)
            for j = 1:Ncls(iC)
                TempConf(i, j) = nanmean(mP_clust(sinds == i, sinds == j), 'all');
            end
        end
        
        ConfMat{iC} = TempConf;
        ConfOut{iC} = TempConf;
    end
end

function [subsamples, pair_count] = sampleCount(population_size, sample_size, num_samples)
%  random subsamples and count pair occurrences
    subsamples = zeros(num_samples, sample_size);
    for i = 1:num_samples
        subsamples(i, :) = randsample(population_size, sample_size, false);
    end
    
    % Count how many times each pair appears
    pair_count = zeros(population_size, population_size);
    for i = 1:num_samples
        idx = subsamples(i, :);
        pair_count(idx, idx) = pair_count(idx, idx) + 1;
    end
end