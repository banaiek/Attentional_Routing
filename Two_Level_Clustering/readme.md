For running consensus-based two-level clustering, you need to run main_sampling_clustering().
This function uses synthetic data for the clustering, but you can replace Data_Input with your input, which in our case is HFAb_triggered_HFA.
The function first calculates PCA, extracts PC1, and then runs two-level consensus-based clustering on it, and saves the results.
The results include a confusion matrix and a pairwise grouping likelihood matrix for each level, as well as cluster labels for each electrode and each cluster number.
