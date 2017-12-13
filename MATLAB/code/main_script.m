clear;
clc;
close all;

%% Plot I-V curve fits
Na_IV_curve;
K1_IV_curve;
K_IV_curve;
Kp_IV_curve;
LCC_IV_curve;

%% Plot gating parameter fits
Na_gating_m;
Na_gating_h;
Na_gating_j;
K1_gating;
K_gating_X;
K_gating_Xi;
Kp_gating;
LCC_gating_d;

%% Plot f-gate fit
LCC_f_gate;

%% Sodium-calcium exchanger model
% Fit parameters
fit_NCX_parameters;
% Run a simulation of the model for an AP input
NCX_sim;

%% Calculate BG parameters
cardiac_AP_bg_parameters;

%% Find the conserved moieties of the model
cardiac_AP_conservation;

%% Run a simulation of a single cardiac action potential
cardiac_AP_sim_single;

%% Simulate and plot effects of conserved moieties on drift (Figure 4)
% Run simulations of variants under conservative and nonconservative
% stimuli
cardiac_AP_sim_conservative_var;
cardiac_AP_sim_nonconservative_var;

cardiac_AP_sim_conservative;
cardiac_AP_sim_nonconservative;

cardiac_AP_sim_ie_stat_con;
cardiac_AP_sim_ie_stat_ncs;

% Plot a comparison of variants with different stimulus currents
stim_current_comparison;

%% Simulate and plot effects of initial conditions on drift (Figure 5)
% Run simulations of variants with different initial conditions
cardiac_AP_sim_ic2;
cardiac_AP_sim_ic3;

cardiac_AP_sim_ie_ic2;
cardiac_AP_sim_ie_ic3;

cardiac_AP_sim_ic2_var;
cardiac_AP_sim_ic3_var;

% Plot a comparison of variants with different initial conditions
ss_comparison;

%% Calculate the conserved moieties for each variant under the initial conditions
ic_conserved_moieties;

