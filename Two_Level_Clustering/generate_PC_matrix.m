function [Clust_vec, Loading_Vec] = generate_PC_matrix(DataInput, Num_Channels, Num_Timepoints)
% get PC1 coefficient matrix from timeseries data
% Inputs:
%   DataInput - data matrix (channels x channels x time)
%   Num_Channels - number of channels
%   Num_Timepoints - number of time points
% Outputs:
%   Clust_vec - z-scored signed clustering matrix
%   Loading_Vec - PC1 loading vector

    % Reshape data for PCA (all channel pairs x time)
    vecSize = Num_Channels * Num_Channels;
    reshaped_Vec = reshape(DataInput, [vecSize, Num_Timepoints]);
    
    % Perform PCA on all data
    [coeff, ~, ~, ~, explained] = pca(reshaped_Vec');
    fprintf('  PC1 explains %.2f%% of variance\n', explained(1));
    
    % Extract and normalize PC1 loadings
    Loading_Vec = coeff(:, 1);
    Loading_Vec = Loading_Vec / max(abs(Loading_Vec));
    
    % Create clustering matrix (z-scored absolute values with sign preserved)
    rLoading_Vec = reshape(Loading_Vec, [Num_Channels, Num_Channels]);
    Clust_vec = zscore(abs(rLoading_Vec)) .* sign(rLoading_Vec);
end