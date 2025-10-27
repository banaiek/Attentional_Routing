
# Attentional_Routing

These scripts are related to the paper “Fast Attentional Information Routing via High-Frequency Bursts in the Human Brain” by Banaie Boroujeni et al.

Written by Kianoush Banaie Boroujeni — 2024

The documentation is organized into three main parts:
	1.	Modeling: Implements the four-network spiking neural network model for the target-detection task.
	2.	Burst Detection: Detects LFP bursts.
	3.	Clustering: Computes the temporal correlation of HFAbs between electrodes. For different cluster numbers, it identifies the most reliable clusters using a two-level consensus-based clustering approach, providing outputs such as the confusion matrix and cluster reliability.
	4. A plotting folder which uses source data to plot the results
