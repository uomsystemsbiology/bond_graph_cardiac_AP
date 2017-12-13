clear;
clc;
close all;

%% Options
run_simulation = true;
sim_length = 30*60; % unit seconds
segment_length = 60; % unit seconds

%% Set directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];
sim_result_dir = [storage_dir 'conservative_stim_result' filesep];
if ~exist(sim_result_dir,'dir')
    mkdir(sim_result_dir);
end

%% Run simulation
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
            [VOI, STATES] = cardiac_AP_exported(tspan,'con');
        else
            [VOI, STATES] = cardiac_AP_exported(tspan,'con',X0);
        end
        toc
        fprintf(['Simulated segment ' num2str(i_segment) ' of ' num2str(num_segments) '\n']);
        
        X0 = STATES(end,:);
        if i_segment > 1
            STATES = single(STATES);
        end
        save([sim_result_dir 'AP_sim_states_con_' num2str(i_segment) '.mat'],'VOI','STATES');   
    end
end

