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

cNao_st = 140; % unit mM
cNai_st = 10; % unit mM
cKo_st = 5.4;
cKi_st = 145;

N_A = 6.022e23;

A_cap = 1.534e-4; % Unit cm^2

x_LCC = 50000/N_A*1e15; % unit fmol
x_Na_channel = 122720/N_A*1e15; % unit fmol
x_K1_channel = 4261/N_A*1e15; % unit fmol
x_K_channel = 5369/N_A*1e15; % unit fmol
x_Kp_channel = 725/N_A*1e15; % unit fmol

%% Define volumes (unit pL)
W_i = 38;
W_e = 5.182;

%% Import IV parameters
% Na channel
G_Na = 16*A_cap; % Unit mA/V
E_Na = R*T/F*log(cNao_st/cNai_st);
E_Na_norm = F*E_Na/(R*T);
G_GHK = 2*G_Na*(1-exp(E_Na_norm))/(cNai_st - cNao_st*exp(E_Na_norm))*R*T/F; % Unit mA/mM
P_Na = G_GHK/F * 1e12; % Unit pL/s

% L-type calcium channel
% Load P_GHK (unit pL/s)
load([storage_dir 'LCC_P_GHK_Ca.mat']);
P_LCC_Ca = P_GHK; % Unit pL/s
load([storage_dir 'LCC_P_GHK_K.mat']);
P_LCC_K = P_GHK; % Unit pL/s

% K1 channel
load([storage_dir 'K1_G_GHK.mat']);
P_K1 = G_GHK/F * 1e12; % Unit pL/s

% K channel
% Load G_GHK
load([storage_dir 'K_G_GHK.mat']);
P_K = G_GHK/F * 1e12; % Unit pL/s

% Kp channel
load([storage_dir 'Kp_G_GHK.mat']);
P_Kp = G_GHK/F * 1e12; % Unit pL/s

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

%% Load Na/K pump kinetic parameters
load([data_dir 'NaK_kinetic_fitting_results.mat']);
params = params_vec;

fast_kinetic_constant = 1e6;

K_d_Nai = params(3);
k_3b_p = fast_kinetic_constant;
k_3b_m = fast_kinetic_constant*(0.5*K_d_Nai);

k_4b_p = fast_kinetic_constant;
k_4b_m = fast_kinetic_constant*(2*K_d_Nai);

K_d_Nai_0 = params(1);
k_5b_m_0 = fast_kinetic_constant;
k_5b_p_0 = fast_kinetic_constant/K_d_Nai_0;

k_6b_p = params(7);
k_6b_m = k_6b_p/6.3;

k_7b_p = params(8);
k_7b_m = params(11);

K_d_Nae_0 = params(2);
k_8b_m_0 = fast_kinetic_constant;
k_8b_p_0 = fast_kinetic_constant*K_d_Nae_0;

K_d_Nae = params(4);
k_9b_m = fast_kinetic_constant;
k_9b_p = fast_kinetic_constant*2*K_d_Nae;

k_10b_m = fast_kinetic_constant;
k_10b_p = fast_kinetic_constant*0.5*K_d_Nae;

K_d_Ke = params(5);
k_11b_p = fast_kinetic_constant;
k_11b_m = fast_kinetic_constant*0.5*K_d_Ke;

k_12b_p = fast_kinetic_constant;
k_12b_m = fast_kinetic_constant*2*K_d_Ke;

k_13b_p = params(9);
k_13b_m = params(12);

K_d_MgATP = params(6);
k_14b_p = fast_kinetic_constant;
k_14b_m = fast_kinetic_constant*K_d_MgATP;

k_15b_p = params(10);
k_15b_m = params(13);

% Calculate remaining parameter
G_MgATP_0 = 11900;
RT = 8.314*310;
K_MgATP = exp(-G_MgATP_0/RT)*10^6;
K_d_Ki = sqrt((k_3b_m*k_4b_m*k_5b_m_0*k_6b_m*k_7b_m*k_8b_m_0*k_9b_m*k_10b_m*k_11b_m*k_12b_m*k_13b_m*k_14b_m*k_15b_m*K_MgATP)/(k_3b_p*k_4b_p*k_5b_p_0*k_6b_p*k_7b_p*k_8b_p_0*k_9b_p*k_10b_p*k_11b_p*k_12b_p*k_13b_p*k_14b_p*k_15b_p));

k_1b_m = fast_kinetic_constant;
k_1b_p = fast_kinetic_constant*2*K_d_Ki;

k_2b_m = fast_kinetic_constant;
k_2b_p = fast_kinetic_constant*0.5*K_d_Ki;

Delta_NaK = params(14);

%% Calculate NCX kinetic parameters (in terms of concentration)
load([storage_dir 'NCX_fitting_results.mat']);
kappa_fast = 1e6;
kappa_3 = params_vec(1);
kappa_6 = params_vec(2);
K_1 = params_vec(3);
K_2 = params_vec(4);
K_3 = params_vec(5);
K_4 = params_vec(6);
K_5 = params_vec(7);
K_6 = params_vec(8);

K_Nai = W_i*params_vec(9);
K_Cai = W_i*params_vec(10);
K_Nae = K_Nai;
K_Cae = K_Cai;

Delta_NCX = params_vec(11);

kf_NCX = [kappa_fast*K_1;...
    kappa_fast*K_2*K_Cai;...
    kappa_3*K_3;...
    kappa_fast*K_4;...
    kappa_fast*K_5*K_Nae^3;...
    kappa_6*K_6];
kr_NCX = [kappa_fast*K_2*K_Nai^3;...
    kappa_fast*K_3;...
    kappa_3*K_4;...
    kappa_fast*K_5*K_Cae;...
    kappa_fast*K_6;...
    kappa_6*K_1];

%% Kinetics of calcium buffering
fast_kinetic_constant = 1e6;

Km_TRPN = 0.5e-3; % unit mM
Km_CMDN = 2.38e-3; % unit mM

kf_TRPN = fast_kinetic_constant;
kr_TRPN = fast_kinetic_constant*Km_TRPN;

kf_CMDN = fast_kinetic_constant;
kr_CMDN = fast_kinetic_constant*Km_CMDN;

kf_Ca_buffer = [kf_TRPN; kf_CMDN];
kr_Ca_buffer = [kr_TRPN; kr_CMDN];

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
num_rows = 73;

N_f = [];
N_r = [];

for i_system = 1:num_subsystems
    T = calcT(struct_stoich{i_system}.I_vec,num_rows);

    struct_stoich{i_system}.T = T;
    
    N_f = [N_f T*struct_stoich{i_system}.N_f];
    N_r = [N_r T*struct_stoich{i_system}.N_r];
end

N_fT = transpose(N_f);
N_rT = transpose(N_r);

N = N_r - N_f;
N_T = N_rT - N_fT;

num_cols = size(N,2);
I = eye(num_cols);

M = [I N_fT; I N_rT];
M_rref = rref(M);

M_T = transpose(M);
M_T_rref = rref(M_T);

%% Set up the vectors

% K1 channel
alpha_K1_bg = params_K1(1)*1e3; % unit s^-1
beta_K1_bg = params_K1(3)*1e3; % unit s^-1

kf_K1 = [P_K1/x_K1_channel/sqrt(cKo_st); ... % R_GHK
    alpha_K1_bg]; % Rg

kr_K1 = [P_K1/x_K1_channel/sqrt(cKo_st); ... % R_GHK
    beta_K1_bg]; % Rg

% K channel
alpha_X_bg = params_X(1); % unit s^-1
beta_X_bg = params_X(3); % unit s^-1

alpha_Xi_bg = params_Xi(1)*1e3; % unit s^-1
beta_Xi_bg = params_Xi(3)*1e3; % unit s^-1

kf_K = [P_K/x_K_channel/sqrt(cKo_st); ... % R_GHK
    2*alpha_X_bg; ... % RX00
    alpha_X_bg; ... % RX00
    2*alpha_X_bg; ... % RX00
    alpha_X_bg; ... % RX00
    alpha_Xi_bg; ... % RX00
    alpha_Xi_bg; ... % RX00
    alpha_Xi_bg]; % RX00

kr_K = [P_K/x_K_channel/sqrt(cKo_st); ... % R_GHK
    beta_X_bg; ... % RX00
    2*beta_X_bg; ... % RX00
    beta_X_bg; ... % RX00
    2*beta_X_bg; ... % RX00
    beta_Xi_bg; ... % RX00
    beta_Xi_bg; ... % RX00
    beta_Xi_bg]; % RX00

% Kp channel
alpha_Kp_bg = params_Kp(1)*1e3; % unit s^-1
beta_Kp_bg = params_Kp(3)*1e3; % unit s^-1

kf_Kp = [P_Kp/x_Kp_channel; ... % R_GHK
    alpha_Kp_bg]; % Rg

kr_Kp = [P_Kp/x_Kp_channel; ... % R_GHK
    beta_Kp_bg]; % Rg

% Na channel
alpha_m0_bg = params_m(1)*1e3; % unit s^-1
beta_m0_bg = params_m(3)*1e3; % unit s^-1

alpha_h0_bg = params_h(1)*1e3; % unit s^-1
beta_h0_bg = params_h(3)*1e3; % unit s^-1

alpha_j0_bg = params_j(1)*1e3; % unit s^-1
beta_j0_bg = params_j(3)*1e3; % unit s^-1

kf_Na = [P_Na/x_Na_channel; ... % R_GHK
    3*alpha_m0_bg; ... % Rm000
    3*alpha_m0_bg; ... % Rm001
    3*alpha_m0_bg; ... % Rm010
    3*alpha_m0_bg; ... % Rm011
    2*alpha_m0_bg; ... % Rm100
    2*alpha_m0_bg; ... % Rm101
    2*alpha_m0_bg; ... % Rm110
    2*alpha_m0_bg; ... % Rm111
    alpha_m0_bg; ... % Rm200
    alpha_m0_bg; ... % Rm201
    alpha_m0_bg; ... % Rm210
    alpha_m0_bg; ... % Rm211
    alpha_h0_bg; ... % Rh000
    alpha_h0_bg; ... % Rh001
    alpha_h0_bg; ... % Rh100
    alpha_h0_bg; ... % Rh101
    alpha_h0_bg; ... % Rh200
    alpha_h0_bg; ... % Rh201
    alpha_h0_bg; ... % Rh300
    alpha_h0_bg; ... % Rh301
    alpha_j0_bg; ... % Rj000
    alpha_j0_bg; ... % Rj010
    alpha_j0_bg; ... % Rj100
    alpha_j0_bg; ... % Rj110
    alpha_j0_bg; ... % Rj200
    alpha_j0_bg; ... % Rj210
    alpha_j0_bg; ... % Rj300
    alpha_j0_bg]; % Rj310

kr_Na = [P_Na/x_Na_channel; ... % R_GHK
    beta_m0_bg; ... % Rm000
    beta_m0_bg; ... % Rm001
    beta_m0_bg; ... % Rm010
    beta_m0_bg; ... % Rm011
    2*beta_m0_bg; ... % Rm100
    2*beta_m0_bg; ... % Rm101
    2*beta_m0_bg; ... % Rm110
    2*beta_m0_bg; ... % Rm111
    3*beta_m0_bg; ... % Rm200
    3*beta_m0_bg; ... % Rm201
    3*beta_m0_bg; ... % Rm210
    3*beta_m0_bg; ... % Rm211
    beta_h0_bg; ... % Rh000
    beta_h0_bg; ... % Rh001
    beta_h0_bg; ... % Rh100
    beta_h0_bg; ... % Rh101
    beta_h0_bg; ... % Rh200
    beta_h0_bg; ... % Rh201
    beta_h0_bg; ... % Rh300
    beta_h0_bg; ... % Rh301
    beta_j0_bg; ... % Rj000
    beta_j0_bg; ... % Rj010
    beta_j0_bg; ... % Rj100
    beta_j0_bg; ... % Rj110
    beta_j0_bg; ... % Rj200
    beta_j0_bg; ... % Rj210
    beta_j0_bg; ... % Rj300
    beta_j0_bg]; % Rj310

% L-type calcium channel
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

kf_LCC = [P_LCC_Ca/x_LCC; ... % R_GHK_Ca1
    P_LCC_Ca/x_LCC; ... % R_GHK_Ca2
    P_LCC_K/x_LCC; ... % R_GHK_K1
    P_LCC_K/x_LCC; ... % R_GHK_K2
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

kr_LCC = [P_LCC_Ca/x_LCC; ... % R_GHK_Ca1
    P_LCC_Ca/x_LCC; ... % R_GHK_Ca2
    P_LCC_K/x_LCC; ... % R_GHK_K1
    P_LCC_K/x_LCC; ... % R_GHK_K2
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

% Na/K pump
kf_NaK = transpose([k_1b_p k_2b_p k_3b_p k_4b_p k_5b_p_0 k_6b_p k_7b_p k_8b_p_0 ...
    k_9b_p k_10b_p k_11b_p k_12b_p k_13b_p k_14b_p k_15b_p]);
kr_NaK = transpose([k_1b_m k_2b_m k_3b_m k_4b_m k_5b_m_0 k_6b_m k_7b_m k_8b_m_0 ...
    k_9b_m k_10b_m k_11b_m k_12b_m k_13b_m k_14b_m k_15b_m]);

% Overall
kf = [kf_K1; kf_K; kf_Kp; kf_Na; kf_LCC; kf_NaK; kf_NCX; kf_Ca_buffer];
kr = [kr_K1; kr_K; kr_Kp; kr_Na; kr_LCC; kr_NaK; kr_NCX; kr_Ca_buffer];

k = [kf;kr];
W = [ones(size(N,2),1);W_i;W_e;ones(10,1);W_i;W_e;ones(16,1);W_i;W_e;ones(12,1);...
    ones(15,1);W_i;W_i;W_i;W_i;ones(6,1);W_i;W_i;W_i;W_i];

lambdaW = exp(pinv(M)*log(k));
lambda = lambdaW./W;
kappa = lambda(1:size(N,2));
K = lambda(size(N,2)+1:end);

save([output_dir 'cardiac_AP_parameters.mat'],'params_m','params_h','params_j','params_d','params_f',...
    'params_K1','params_X','params_Xi','params_Kp','kappa','K','Delta_NaK','Delta_NCX');

%% Checks
N_rref = rref(N);
R_mat = null(N,'r');

K_eq = kf./kr;
zero_est = R_mat'*K_eq;

k_est = exp(M*log(lambdaW));
diff = sum(abs((k-k_est)./k));
