clear;
clc;
close all;

%% Options
run_optimisation = false;

%% Set directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];

%% Steady-state gating parameters and time constants

V = transpose(-120:1:60); % unit mV

alpha_m = 0.31*(V+47.13)./(1-exp(-0.1*(V+47.13)));
beta_m = 0.08*exp(-V/11);

tau_m = 1./(alpha_m + beta_m);
g_ss_m = alpha_m./(alpha_m + beta_m);


%% Fit bond graph parameters to model
% params: [kf, zf, kr, zr];
error_func_alpha = @(params) square_error(alpha_m - calc_alpha(params,V/1000));
error_func_beta = @(params) square_error(beta_m - calc_beta(params,V/1000));
error_func_gss = @(params) square_error(g_ss_m - p2gss(params,V));
error_func_tau = @(params) square_error(tau_m- p2tau(params,V));

error_func = @(params) error_func_alpha(params) + error_func_beta(params) + error_func_gss(params) + error_func_tau(params);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [0; -10; 0; -10];
ub = [Inf; 10; Inf; 10];

options_unc = optimoptions('fminunc','MaxFunEvals',10000);
options_ps = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon);

if run_optimisation
    [params_vec,fval,exitflag,output] = particleswarm(error_func,4,lb,ub,options_ps);
    save([storage_dir 'Na_m_parameters.mat'],'params_vec');
else
    load([storage_dir 'Na_m_parameters.mat']);
end
% [params_vec,fval,exitflag,output,grad,hessian] = fminunc(error_func,params_vec,options_unc);

g_ss_fit = p2gss(params_vec,V);
h1 = figure;
plot(V,g_ss_m,'k--','LineWidth',4);
hold on;
plot(V,g_ss_fit,'k','LineWidth',4);
xlabel('Voltage (mV)');
ylabel('m_{ss}');
set(gca,'FontSize',28);
xlim([-120 60]);
set(gca,'XTick',-120:30:60);
set(gca,'YTick',0:0.2:1);
xticklabels({-120,'',-60,'',0,'',60});
yticklabels({0,'','','','',1});
set(gca,'LineWidth',3);
grid on;

tau_fit = p2tau(params_vec,V);
h2 = figure;
plot(V,tau_m,'k--','LineWidth',4);
hold on;
plot(V,tau_fit,'k','LineWidth',4);
% legend('Luo and Rudy','Fitted');
xlabel('Voltage (mV)');
ylabel('\tau_m (ms)');
set(gca,'FontSize',28);
xlim([-120 60]);
set(gca,'XTick',-120:30:60);
% set(gca,'YTick',0:0.2:1);
xticklabels({-120,'',-60,'',0,'',60});
% yticklabels({0,'','','','',1});
set(gca,'LineWidth',3);
set(gca,'xgrid','on');

print_figure(h1,output_dir,'g_ss_m');
print_figure(h2,output_dir,'tau_m');


