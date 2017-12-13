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

%% Define constants
R = 8.314; % unit J/mol/K
T = 310;
F = 96485;

%% Import gate transition parameters
load([storage_dir 'Na_m_parameters.mat']);
params_m = params_vec;

load([storage_dir 'Na_h_parameters.mat']);
params_h = params_vec;

load([storage_dir 'Na_j_parameters.mat']);
params_j = params_vec;

load([storage_dir 'K1_parameters.mat']);
params_K1 = params_vec;

load([storage_dir 'K_X_parameters.mat']);
params_X = params_vec;

% Calculate parameters for X_i gate
A = exp(-56.26/32.1);
K = 1/(1+A*exp(-120/32.1));

alpha0 = K; % Unit ms^-1
beta0 = A*K; % Unit ms^-1

zf = 0;
zr = R*T/32.1/F;
params_Xi = [alpha0; zf*1e3; beta0; zr*1e3];

A = exp(7.488/5.98);
K = 1/(1+A*exp(-60/5.98));

alpha0 = K; % Unit ms^-1
beta0 = A*K; % Unit ms^-1

zf = 0;
zr = -R*T/5.98/F;

params_Kp = [alpha0; zf*1e3; beta0; zr*1e3];

load([storage_dir 'LCC_d_parameters.mat']);
params_d = params_vec;

load([storage_dir 'LCC_f_parameters.mat']);
params_f = params_vec;

zfK1 = params_K1(2);
zrK1 = params_K1(4);

zfX = params_X(2);
zrX = params_X(4);
zfXi = params_Xi(2);
zrXi = params_Xi(4);

zfKp = params_Kp(2);
zrKp = params_Kp(4);

zfm = params_m(2);
zrm = params_m(4);
zfh = params_h(2);
zrh = params_h(4);
zfj = params_j(2);
zrj = params_j(4);

zfd = params_d(2);
zrd = params_d(4);

zff1 = params_f(2);
zrf1 = params_f(4);
zff2 = params_f(6);
zrf2 = params_f(8);
zff3 = zrf1+zff2-zff1-zrf2;

zrCa = 2;

load([data_dir 'NaK_kinetic_fitting_results.mat']);
params = params_vec;
Delta_NaK = params(14);

%% Load stoichiometric matrices
array_names = {'K1';'K';'Kp';'Na';'LCC';'NaK';'NCX';'Ca_buffer'};
num_subsystems = length(array_names);
struct_stoich = cell(num_subsystems,1);
for i_system = 1:num_subsystems
    sys_name = array_names{i_system};
    struct_stoich{i_system}.name = sys_name;
    [N_f,N_r] = load_stoichiometry(data_dir,sys_name);
    struct_stoich{i_system}.N_f = N_f;
    struct_stoich{i_system}.N_r = N_r;
end

struct_stoich{1}.I_vec = 1:4;
struct_stoich{2}.I_vec = [1:2 5:10];
struct_stoich{3}.I_vec = [1:2 11:12];
struct_stoich{4}.I_vec = 13:30;
struct_stoich{5}.I_vec = [31:32 1:2 33:44];
struct_stoich{6}.I_vec = [45:59 1:2 13:14 60:63];
struct_stoich{7}.I_vec = [64:69 13 31 14 32];
struct_stoich{8}.I_vec = [31 70:73];

num_species = 73;

N_f = [];
N_r = [];

for i_system = 1:num_subsystems
    T = calcT(struct_stoich{i_system}.I_vec,num_species);

    struct_stoich{i_system}.T = T;
    
    N_f = [N_f T*struct_stoich{i_system}.N_f];
    N_r = [N_r T*struct_stoich{i_system}.N_r];
end

N = N_r - N_f;
num_reactions = size(N,2);

% Add charge to stoichiometry
N = [N; zeros(1,num_reactions)];
% K1
N(end,1) = -1;
N(end,2) = zrK1 - zfK1;
% K
N(end,3) = -1;
N(end,4:7) = zrX - zfX;
N(end,8:10) = zrXi - zfXi;
% Kp
N(end,11) = -1;
N(end,12) = zrKp - zfKp;
% Na
N(end,13) = -1;
N(end,14:25) = zrm - zfm;
N(end,26:33) = zrh - zfh;
N(end,34:41) = zrj - zfj;
% LCC
N(end,42:43) = -2;
N(end,44:45) = -1;
N(end,46:51) = zrd - zfd;
N(end,52:55) = zrf1 - zff1;
N(end,56:59) = zrf2 - zff2;
N(end,60:63) = -zff3;
% skip 6 - 69
% Na/K
% skip 4 - 73
N(end,74) = Delta_NaK;
% skip 2 - 76
N(end,77) = -(1+Delta_NaK);
% skip 7 - 84
% NCX
% skip 5 - 89
N(end,90) = 1;

%% Calculate conserved moieties
G = transpose(null(transpose(N),'r'));

% Hold ATP, ADP, H and Pi constant
idx_var = setdiff(1:size(N,1),60:63);
G_var_ions = transpose(null(transpose(N(idx_var,:)),'r'));
G_var_ions_chem = sp_consmoiety(N(idx_var,:));
G_var_ions = [G_var_ions_chem; G_var_ions(end,:)];
G_var_ions(end,:) = G_var_ions(end,:) -2*G_var_ions(7,:) +4*G_var_ions(8,:) + 2*G_var_ions(11,:) + 2*G_var_ions(12,:);

% Hold extracellular ion concentrations constant
idx_var = setdiff(1:size(N,1),[2 14 32 60:63]);
G_var_int = transpose(null(transpose(N(idx_var,:)),'r'));

% Hold intracellular and extracellular ion concentrations constant
idx_var = setdiff(1:size(N,1),[1:2 13:14 31:32 60:63]);
G_const_ions = transpose(null(transpose(N(idx_var,:)),'r'));

G_chem = transpose(null(transpose(N(1:end-1,:)),'r'));

save([output_dir 'cardiac_AP_cm.mat'],'G_var_ions','G_var_int','G_const_ions');

