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

%% Plot I-V curves
cKo_st = 5.4;
cKi_st = 145;

A_geo = 1.534e-4; % Unit cm^2
G_Kp = 0.0183*A_geo; % Unit mA/V

V = (-120:1:60)/1000;
E_K = R*T/F*log(cKo_st/cKi_st);
I_lin = G_Kp*(V-E_K);

error_func = @(G_GHK) square_error(I_lin(V>=0) - calc_IGHK(G_GHK,V(V>=0),cKi_st,cKo_st));

A = [];
b = [];
Aeq = [];
beq = [];
lb = [-Inf];
ub = [Inf];

options_ps = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fminunc,'SwarmSize',200);

if run_optimisation
    [G_GHK,fval,exitflag,output] = particleswarm(error_func,1,lb,ub,options_ps);
    save([storage_dir 'Kp_G_GHK.mat'],'G_GHK');
else
    load([storage_dir 'Kp_G_GHK.mat']);
end
I_GHK = calc_IGHK(G_GHK,V,cKi_st,cKo_st);

h = figure;
plot(1000*V,1e6*I_lin,'k--',1000*V,1e6*I_GHK,'k','LineWidth',2);
legend('LRd','BG','Location','southeast');
ylabel('Current (nA)');
xlabel('Voltage (mV)');
set(gca,'FontSize',16);

print_figure(h,output_dir,'Kp_IV_curve');