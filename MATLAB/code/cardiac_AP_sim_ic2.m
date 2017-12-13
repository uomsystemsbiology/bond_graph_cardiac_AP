clear;
clc;
close all;

%% Options
run_simulation = true;
save_figures = false;
sim_length = 10*60; % unit seconds
segment_length = 60; % unit seconds

%% Set directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];
sim_result_dir = [storage_dir 'ic2_e_const' filesep];
if ~exist(sim_result_dir,'dir')
    mkdir(sim_result_dir);
end

%% Run simulation
W_i = 38; % unit pL

% Create time intervals for simulation
t_vec = 0:segment_length:sim_length;
if t_vec(end) ~= sim_length
    t_vec = [t_vec sim_length];
end
num_segments = length(t_vec) - 1;

if run_simulation
    tic
    for i_segment = 1:num_segments
        tspan = t_vec(i_segment):0.001:t_vec(i_segment+1);
        
        if i_segment == 1
            [~, X0] = cardiac_AP_exported(0,'init');
            X0(4) = X0(4) - 1*W_i; % Ki
            X0(6) = X0(6) + 1*W_i; % Ki
            [VOI, STATES] = cardiac_AP_exported(tspan,'con',X0);
        else
            [VOI, STATES] = cardiac_AP_exported(tspan,'con',X0);
        end
        toc
        fprintf(['Simulated segment ' num2str(i_segment) ' of ' num2str(num_segments) '\n']);
        
        X0 = STATES(end,:);
        STATES = single(STATES);
        save([sim_result_dir 'AP_sim_states_' num2str(i_segment) '.mat'],'VOI','STATES');   
    end
end
