clear;
clc;
close all;

%% Options
run_optimisation = false;

%% Set up directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];

%% Define constants
R = 8.314;
T = 310;
F = 96485;

P_NaK = 0.01833;

%% Plot I-V curves
cKo_st = 5.4;
cKo = 5.4;
cKi_st = 145;
cNao_st = 140;
cNai_st = 10;

A_geo = 1.534e-4; % Unit cm^2
G_K = 0.282*A_geo*sqrt(cKo/5.4); % Unit mA/V

V = (-120:1:60)/1000;
E_K = R*T/F*log(cKo/cKi_st);
E_K_LR = R*T/F*log((cKo + P_NaK*cNao_st)/(cKi_st + P_NaK*cNai_st));
I_lin = G_K*(V-E_K_LR); % Unit mA

E_K_norm = F*E_K/(R*T);
G_GHK = 5e-8; % Unit mA/mM
V_norm = F*V/R/T;
GHK_factor = V_norm./(1-exp(-V_norm));
idx_nan = isnan(GHK_factor);
GHK_factor(idx_nan) = 1;
I_GHK = G_GHK*GHK_factor.*(cKi_st - cKo*exp(-V_norm));
% I_GHK = G_GHK*GHK_factor.*(100 - cNao_st*exp(-V_norm));

error_func = @(G_GHK) square_error(I_lin(V>=-0.02 & V<0.03) - calc_IGHK(G_GHK,V(V>=-0.02 & V<0.03),cKi_st,cKo));

A = [];
b = [];
Aeq = [];
beq = [];
lb = [-Inf];
ub = [Inf];


options_ps = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fminunc,'SwarmSize',200);

if run_optimisation
    [G_GHK,fval,exitflag,output] = particleswarm(error_func,1,lb,ub,options_ps);
    save([storage_dir 'K_G_GHK.mat'],'G_GHK');
else
    load([storage_dir 'K_G_GHK.mat']);
end
I_GHK = calc_IGHK(G_GHK,V,cKi_st,cKo);

h = figure;
plot(1000*V,1e6*I_lin,'k--',1000*V,1e6*I_GHK,'k','LineWidth',2);
legend('LRd','BG','Location','southeast');
ylabel('Current (nA)');
xlabel('Voltage (mV)');
set(gca,'FontSize',16);

print_figure(h,output_dir,'K_IV_curve');