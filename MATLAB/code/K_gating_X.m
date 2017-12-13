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

V = transpose(-120:1:60);
alpha_X = 7.19e-2 * (V+30) ./ (1-exp(-0.148*(V+30))); % unit s^-1
idx_nan = isnan(alpha_X);
alpha_X(idx_nan) = 7.19e-2 / 0.148;
beta_X = 1.31e-1 * (V+30) ./ (-1+exp(0.0687*(V+30))); % unit s^-1
idx_nan = isnan(beta_X);
beta_X(idx_nan) = 1.31e-1 / 0.0687;

tau_X = 1./(alpha_X + beta_X);
g_ss_X = alpha_X./(alpha_X + beta_X);


%% Fit bond graph parameters to model
% params: [kf, zf, kr, zr];
multiplier = 1 + 4*(V >= 0);
error_func_alpha = @(params) square_error(multiplier.*(alpha_X - calc_alpha(params,V/1000)));
error_func_beta = @(params) square_error(multiplier.*(beta_X - calc_beta(params,V/1000)));
error_func_gss = @(params) square_error(multiplier.*(g_ss_X - p2gss(params,V)));
error_func_tau = @(params) square_error(multiplier.*(tau_X - p2tau(params,V)));

error_func = @(params) error_func_alpha(params) + error_func_beta(params) + 100*error_func_gss(params) + error_func_tau(params);
% error_func = @(params) error_func_alpha(params) + error_func_beta(params) + error_func_gss(params) + error_func_tau(params);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [0; -10; 0; -10];
ub = [Inf; 10; Inf; 10];

options_unc = optimoptions('fminunc','MaxFunEvals',10000);
options_ps = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon,'SwarmSize',200);

if run_optimisation
    [params_vec,fval,exitflag,output] = particleswarm(error_func,4,lb,ub,options_ps);
    save([storage_dir 'K_X_parameters.mat'],'params_vec');
else
    load([storage_dir 'K_X_parameters.mat']);
end

g_ss_fit = p2gss(params_vec,V);
h1 = figure;
plot(V,g_ss_X,'k--','LineWidth',4);
hold on;
plot(V,g_ss_fit,'k','LineWidth',4);
% legend('Luo and Rudy','Fitted','Location','northwest');
xlabel('Voltage (mV)');
ylabel('X_{ss}');
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
plot(V,1000*tau_X,'k--','LineWidth',4);
hold on;
plot(V,1000*tau_fit,'k','LineWidth',4);
% legend('Luo and Rudy','Fitted','Location','northwest');
xlabel('Voltage (mV)');
ylabel('\tau_X (ms)');
set(gca,'FontSize',28);
xlim([-120 60]);
set(gca,'XTick',-120:30:60);
% set(gca,'YTick',0:0.2:1);
xticklabels({-120,'',-60,'',0,'',60});
% yticklabels({0,'','','','',1});
set(gca,'LineWidth',3);
set(gca,'xgrid','on');

print_figure(h1,output_dir,'g_ss_X');
print_figure(h2,output_dir,'tau_X');



