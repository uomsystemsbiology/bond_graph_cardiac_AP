clear;
clc;

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

N_A = 6.022e23;

x_LCC = 50000/N_A*1e15; % unit fmol

%% Define volumes (unit pL)
W_i = 38;
W_e = 5.182;

%% Import IV parameters
A_geo = 1.534e-4; % Unit cm^2
% Load P_GHK (unit pL/s)
load([storage_dir 'LCC_P_GHK_Ca.mat']);
P_Ca = P_GHK;

load([storage_dir 'LCC_P_GHK_K.mat']);
P_K = P_GHK;

%% Import gate transition parameters
load([storage_dir 'LCC_d_parameters.mat']);
params_d = params_vec;

load([storage_dir 'LCC_f_parameters.mat']);
params_f = params_vec;

%% Load stoichiometric matrices
stoich_path = [data_dir 'LCC_forward_matrix.txt'];
stoich_file_id = fopen(stoich_path,'r');

stoich_file_data = textscan(stoich_file_id,'%s','delimiter','\n');
fclose(stoich_file_id);

num_rows = length(stoich_file_data{1});
num_cols = sum(stoich_file_data{1}{1} == ',')+1;

N_f = zeros(num_rows,num_cols);

for i_row = 1:num_rows
    line_str = stoich_file_data{1}{i_row};
    line_split = regexp(line_str,',','split');
    
    for i_col = 1:num_cols
        N_f(i_row,i_col) = str2double(line_split{i_col});
    end
end

N_fT = transpose(N_f);

stoich_path = [data_dir 'LCC_reverse_matrix.txt'];
stoich_file_id = fopen(stoich_path,'r');

stoich_file_data = textscan(stoich_file_id,'%s','delimiter','\n');
fclose(stoich_file_id);

num_rows = length(stoich_file_data{1});
num_cols = sum(stoich_file_data{1}{1} == ',')+1;

N_r = zeros(num_rows,num_cols);

for i_row = 1:num_rows
    line_str = stoich_file_data{1}{i_row};
    line_split = regexp(line_str,',','split');
    
    for i_col = 1:num_cols
        N_r(i_row,i_col) = str2double(line_split{i_col});
    end
end

N_rT = transpose(N_r);

N = N_r - N_f;
N_T = N_rT - N_fT;

I = eye(num_cols);

M = [I N_fT; I N_rT];
M_rref = rref(M);

M_T = transpose(M);
M_T_rref = rref(M_T);

%% Set up the vectors
alpha_d0_bg = params_d(1)*1e3; % unit s^-1
beta_d0_bg = params_d(3)*1e3; % unit s^-1

alpha_f1_0_bg = params_f(1)*1e3; % unit s^-1
beta_f1_0_bg = params_f(3)*1e3; % unit s^-1

alpha_f2_0_bg = params_f(5)*1e3; % unit s^-1
beta_f2_0_bg = params_f(7)*1e3; % unit s^-1

K_f3_0 = beta_f1_0_bg*alpha_f2_0_bg/alpha_f1_0_bg/beta_f2_0_bg; % dimensionless
rate_f3 = 1e5; % unit s^-1

KmCa = 0.6e-3; % Unit mM
rate_fCa = 1e5; % unit s^-1

kf = [P_Ca/x_LCC; ... % R_GHK_Ca1
    P_Ca/x_LCC; ... % R_GHK_Ca2
    P_K/x_LCC; ... % R_GHK_K1
    P_K/x_LCC; ... % R_GHK_K2
    alpha_d0_bg; ... % Rd000
    alpha_d0_bg; ... % Rd010
    alpha_d0_bg; ... % Rd020
    alpha_d0_bg; ... % Rd001
    alpha_d0_bg; ... % Rd011
    alpha_d0_bg; ... % Rd021
    alpha_f1_0_bg; ... % Rf1_000
    alpha_f1_0_bg; ... % Rf1_100
    alpha_f1_0_bg; ... % Rf1_001
    alpha_f1_0_bg; ... % Rf1_101
    alpha_f2_0_bg; ... % Rf2_000
    alpha_f2_0_bg; ... % Rf2_100
    alpha_f2_0_bg; ... % Rf2_001
    alpha_f2_0_bg; ... % Rf2_101
    K_f3_0*rate_f3; ... % Rf3_010
    K_f3_0*rate_f3; ... % Rf3_110
    K_f3_0*rate_f3; ... % Rf3_011
    K_f3_0*rate_f3; ... % Rf3_111
    rate_fCa; ... % RfCa000
    rate_fCa; ... % RfCa100
    rate_fCa; ... % RfCa010
    rate_fCa; ... % RfCa110
    rate_fCa; ... % RfCa020
    rate_fCa]; % RfCa120

kr = [P_Ca/x_LCC; ... % R_GHK_Ca1
    P_Ca/x_LCC; ... % R_GHK_Ca2
    P_K/x_LCC; ... % R_GHK_K1
    P_K/x_LCC; ... % R_GHK_K2
    beta_d0_bg; ... % Rd000
    beta_d0_bg; ... % Rd010
    beta_d0_bg; ... % Rd020
    beta_d0_bg; ... % Rd001
    beta_d0_bg; ... % Rd011
    beta_d0_bg; ... % Rd021
    beta_f1_0_bg; ... % Rf1_000
    beta_f1_0_bg; ... % Rf1_100
    beta_f1_0_bg; ... % Rf1_001
    beta_f1_0_bg; ... % Rf1_101
    beta_f2_0_bg; ... % Rf2_000
    beta_f2_0_bg; ... % Rf2_100
    beta_f2_0_bg; ... % Rf2_001
    beta_f2_0_bg; ... % Rf2_101
    rate_f3; ... % Rf3_010
    rate_f3; ... % Rf3_110
    rate_f3; ... % Rf3_011
    rate_f3; ... % Rf3_111
    rate_fCa/KmCa; ... % RfCa000
    rate_fCa/KmCa; ... % RfCa100
    rate_fCa/KmCa; ... % RfCa010
    rate_fCa/KmCa; ... % RfCa110
    rate_fCa/KmCa; ... % RfCa020
    rate_fCa/KmCa]; % RfCa120

k = [kf;kr];
W = [ones(28,1);W_i;W_e;W_i;W_e;ones(12,1)];

lambdaW = exp(pinv(M)*log(k));
lambda = lambdaW./W;
kappa = lambda(1:28);
K = lambda(29:end);

save('LR_LCC_ch_parameters.mat','params_d','params_f',...
    'kappa','K');

%% Checks
N_rref = rref(N);
R_mat = null(N,'r');

K_eq = kf./kr;
zero_est = R_mat'*K_eq;

k_est = exp(M*log(lambdaW));
diff = sum(abs((k-k_est)./k));
