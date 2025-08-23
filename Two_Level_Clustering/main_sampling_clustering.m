function main_sampling_clustering()
%  consensus-based two level sampling clustering 
    clear; close all;
    
    % for illustration it uses synthetic timeseries data (100x100x1000) [ you can replace
    % this with HFAb-triggered-HFA between electrodes (DataInput = data,
    % Num_Channels = number of channels, Num_Timepoints =  number of time
    % points for each HFAb-triggered-HFA).

    Num_Channels = 100;
    Num_Timepoints = 1000;   
    DataInput = randn(Num_Channels, Num_Channels, Num_Timepoints);
    
    % Normalize the data
    DataInput = (DataInput - mean(DataInput, 3)) ./ std(DataInput, [], 3);
    
    % parameters
    params.Ncls = 2:8;                  % number of clusters to test
    params.proportion = 0.25;           % subsampling proportion
    params.num_samples = 1000;          % number of iterations
    params.Corr_type = 'correlation';   % distance metric
    
    % Get PC1 coefficient matrix
    fprintf('Generating PC1 coefficient matrix...\n');
    [Clust_vec, Loading_Vec] = generate_PC_matrix(DataInput, Num_Channels, Num_Timepoints);
    
    % Perform consensus-based two level clustering
    fprintf('Performing two-level sampling-based clustering...\n');
    RESClust = perform_sampling_clustering(Clust_vec, params.proportion, ...
                                          params.num_samples, params.Ncls, ...
                                          params.Corr_type);
    
    RESClust.Clust_vec = Clust_vec;
    RESClust.Loading_Vec = Loading_Vec;
    RESClust.params = params;
    RESClust.Num_Channels = Num_Channels;
    
    save('Clustering_RES_Synthetic.mat', 'RESClust', 'params', '-v7.3');
    fprintf('Results saved to Clustering_RES_Synthetic.mat\n');
    
    fprintf('\nClustering completed!\n');
    fprintf('Tested %d cluster configurations (k=%d to %d)\n', ...
            length(params.Ncls), params.Ncls(1), params.Ncls(end));
    fprintf('Used %d iterations per configuration\n', params.num_samples);
end