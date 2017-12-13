clear;
clc;
close all;

%% Set up directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];

%% Define constants
z = 2;
R = 8.314;
T = 310;
F = 96485;

gamma_Cai = 1;
gamma_Cao = 0.341;
gamma_Ki = 0.75;
gamma_Ko = 0.75;

%% Plot I-V curves
cCao_st = 1.8;
cCai_st = 0.12e-3;

A_cap = 1.534e-4; % Unit cm^2
P_Ca = 5.4e-4*A_cap; % Unit cm^3/s

V = (-120:1:60)/1000;
V_norm = z*F*V/R/T;
I_LR = P_Ca * z*F*V_norm .* (gamma_Cai * cCai_st * exp(V_norm) - gamma_Cao* cCao_st) ./ (exp(V_norm) - 1) ; % Unit uA

idx_nan = find(isnan(I_LR));
I_LR(idx_nan) = P_Ca * z*F * (gamma_Cai * cCai_st * exp(V_norm(idx_nan)) - gamma_Cao* cCao_st);

P_GHK = P_Ca*gamma_Cao; % Unit cm^3/s
I_GHK = P_GHK * z*F*V_norm .* (cCai_st * exp(V_norm) - cCao_st) ./ (exp(V_norm) - 1) ; % Unit uA
idx_nan = find(isnan(I_GHK));
I_GHK(idx_nan) = P_GHK * z*F* (cCai_st * exp(V_norm(idx_nan)) - cCao_st);

h = figure;
plot(1000*V,1000*I_LR,'k--',1000*V,1000*I_GHK,'k','LineWidth',2);
legend('LRd','BG','Location','southeast');
ylabel('Current (nA)');
xlabel('Voltage (mV)');
set(gca,'FontSize',16);

P_GHK = P_GHK*1e9; % Unit pL/s

print_figure(h,output_dir,'LCC_Ca_IV_curve');
save([storage_dir 'LCC_P_GHK_Ca.mat'],'P_GHK');

%% Plot I-V curves (K current)
cKo_st = 5.4;
cKi_st = 145;

P_K = 1.93e-7*A_cap; % Unit cm^3/s

V = (-120:1:60)/1000;
V_norm = z*F*V/R/T;
I_LR = P_K * z*F*V_norm .* (gamma_Ki * cKi_st * exp(V_norm) - gamma_Ko* cKo_st) ./ (exp(V_norm) - 1) ; % Unit uA

idx_nan = find(isnan(I_LR));
I_LR(idx_nan) = P_K * z*F * (gamma_Ki * cKi_st * exp(V_norm(idx_nan)) - gamma_Ko* cKo_st);

P_GHK = P_K*gamma_Ko; % Unit cm^3/s
I_GHK = P_GHK * z*F*V_norm .* (cKi_st * exp(V_norm) - cKo_st) ./ (exp(V_norm) - 1) ; % Unit uA
idx_nan = find(isnan(I_GHK));
I_GHK(idx_nan) = P_GHK * z*F* (cKi_st * exp(V_norm(idx_nan)) - cKo_st);

h = figure;
plot(1000*V,1000*I_LR,'k--',1000*V,1000*I_GHK,'k','LineWidth',2);
legend('LRd','BG','Location','southeast');
ylabel('Current (nA)');
xlabel('Voltage (mV)');
set(gca,'FontSize',16);

P_GHK = P_GHK*1e9; % Unit pL/s

print_figure(h,output_dir,'LCC_K_IV_curve');
save([storage_dir 'LCC_P_GHK_K.mat'],'P_GHK');