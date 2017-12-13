clear;
clc;
close all;

%% Set directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];

%% Define constants and options
R = 8.314;
T = 310;
F = 96485;

print_figures = true;

%% Load stoichiometric information and constants
% load([storage_dir 'fitting_results.mat']);

syms Delta
% Delta_num = params_vec(21);

% Load forward matrix
stoich_path = [data_dir 'NCX_forward.csv'];
N_f = csvread(stoich_path);
num_re = size(N_f,2);
N_f_sym = sym([N_f; zeros(1,num_re)]);
N_f = [N_f; zeros(1,num_re)];

% N_f(end,5) = -Delta_num;
% N_f(end,8) = 1+Delta_num;

N_f_sym(end,6) = Delta;

% Load reverse matrix
stoich_path = [data_dir 'NCX_reverse.csv'];
N_r = csvread(stoich_path);
N_r = [N_r; zeros(1,num_re)];
N_r_sym = sym(N_r);

N_r_sym(end,6) = 1+Delta;

% Calculate stoichiometric matrix
N = N_r - N_f;
N_sym = N_r_sym - N_f_sym;

% Define indices for rapid equilibrium approximation
Idx_fast = [1 2 4 5];
Idx_slow = [3 6];

Idx_e = 7:11;
Idx_n = [1 4];
Idx_d = setdiff(setdiff(1:11,Idx_e),Idx_n);

W_i = 38;
W_e = 5.182;
W_total = W_i + W_e;

%% Generate code for cycling velocity
subscript_labels = {'1';'2';'3';'4';'5';'6';'Nai';'Cai';'Nae';'Cae';'mem'};

reaction_labels = {'1';'2';'3';'4';'5';'6'};

struct_input = struct('N_f',N_f_sym,...
    'N_r',N_r_sym,...
    'N',N_sym,...
    'Idx_e',Idx_e,...
    'Idx_n',Idx_n,...
    'Idx_d',Idx_d,...
    'Idx_fast',Idx_fast,...
    'Idx_slow',Idx_slow,...
    'subscript_labels',{subscript_labels},...
    'reaction_labels',{reaction_labels});

%% Load NCX parameters and simulation conditions
W_i = 38;
W_e = 5.182;
W_total = W_i + W_e;

R = 8.314;
T = 310;
F = 96485;

A_cap = 1.534e-4; % unit cm^2
mem_cap = 1e9 * A_cap; % unit fF

load([storage_dir 'NCX_fitting_results.mat']);

struct_input = struct('Nai',10,...
    'Cae',1.8,...
    'Cai',120e-6,...
    'V',-0.06,...
    'Nae',140);

%% Run simulation of NCX, using action potential waveform
NCX_density = 170; % unit um^-2

% Load AP waveform from Faber and Rudy (2000)
CellML_data = csvread([data_dir 'faber_rudy_2000_data.csv'],1,0);

t_LRd = CellML_data(:,1); % unit ms
V_LRd = CellML_data(:,5); % unit mV
Cai_LRd = CellML_data(:,4); % unit mM

t_start = -100;

% Restrict to an action potential
t_LRd = t_LRd(t_LRd <= 400) + t_start;
V_LRd = V_LRd(t_LRd <= 400);

v_NCX = zeros(length(t_LRd),1);
for i_time = 1:length(t_LRd)
    V_m = V_LRd(i_time)/1000;
    Cai = Cai_LRd(i_time);
    struct_input.V = V_m;
    struct_input.Cai = Cai;
    v_NCX(i_time) = NCX_vss_fitting(params_vec,struct_input,1);
end

I_NCX = -NCX_density*v_NCX*1.6e-5; % unit uA/cm^2

h1 = figure;
plot(t_LRd,I_NCX,'LineWidth',3);
xlabel('Time (ms)');
ylabel('Current density ({\mu}A/cm^2)');
set(gca,'FontSize',18);

print_figure(h1,output_dir,'AP_sim');
