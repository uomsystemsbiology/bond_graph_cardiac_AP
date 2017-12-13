clear;
clc;
close all;

%% Script options
run_optimisation = false;
marker_size = 30;
num_runs = 50;
swarm_size = 1000;
save_figures = true;

%% Set up directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];

%% Load parameters
load([storage_dir 'NCX_fitting_results.mat']);

%% Look at saturation
struct_input = struct('Nai',8,...
    'Cae',1,...
    'Cai',1e-3,...
    'V',0,...
    'Nae',150);

V_vec = (-500:1:500)/1000;

r_I_model = zeros(1,length(V_vec));

for i_V = 1:length(V_vec)
    struct_input.V = V_vec(i_V);
    r_I_model(i_V) = NCX_vss_fitting(params_vec,struct_input,1);
end

h1 = figure;
plot(1000*V_vec,r_I_model,'LineWidth',2);