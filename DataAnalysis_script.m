% This script runs the analyses to reproduce the results and figures from 
% "Motor cortex retains and reorients neural dynamics during motor imagery" 
% by Dekleva et al
% 
% There are some dependencies across scripts, so they may need to be executed
% in the given order. 

%% Load the data 
% i.e. load('ActionImageryData.mat')

%% Figure 2a-e, Figure 3b-d, Supplementary Figure 3
Split_Subspaces_script; 

%% Supplementary Figure 5
Channel_Contributions_script; 

%% Figure 1c, Supplementary Figure 1
Forces_script;

%% Figure 1d
Total_modulations_script;

%% Figure 1f
Single_channel_script;

%% Figure 3a
Monte_carlo_script;

%% Figure 2f,g and Figure 5f,g
Dimensionalities_script;

%% Figure 5c,d
Force_dependence_script;

%% Figure 4
Novel_common_dynamics_script;

%% Figure 6
Decoding_script;

%% Supplementary Figure 6
Visual_control_script;

%% Supplementary Figure 2
EMG_script;
