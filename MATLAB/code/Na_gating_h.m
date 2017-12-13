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
alpha_h = zeros(length(V),1);
beta_h = zeros(length(V),1);

for i_V = 1:length(V)
    if V(i_V) >= -40
        alpha_h(i_V) = 0;
        beta_h(i_V) = 1/(0.13*(1+exp((V(i_V)+10.66)/-11.1)));
    else
        alpha_h(i_V) = 0.135*exp((80+V(i_V))/-6.8);
        beta_h(i_V) = 3.56*exp(0.079*V(i_V)) + 3.1e5 * exp(0.35*V(i_V));
    end
end

tau_h = 1./(alpha_h + beta_h);
g_ss_h = alpha_h./(alpha_h + beta_h);

%% Fit bond graph parameters to model
% params: [kf, zf, kr, zr];
error_func_alpha = @(params) square_error(alpha_h - calc_alpha(params,V/1000));
error_func_beta = @(params) square_error(beta_h - calc_beta(params,V/1000));
error_func_gss = @(params) square_error(g_ss_h - p2gss(params,V));
error_func_tau = @(params) square_error(tau_h- p2tau(params,V));

error_func = @(params) error_func_alpha(params) + error_func_beta(params) + error_func_gss(params) + error_func_tau(params);
% error_func = @(params) error_func_gss(params) + error_func_tau(params);


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
    save([storage_dir 'Na_h_parameters.mat'],params_vec);
else
    load([storage_dir 'Na_h_parameters.mat']);
end
% [params_vec,fval,exitflag,output,grad,hessian] = fminunc(error_func,params_vec,options_unc);

g_ss_fit = p2gss(params_vec,V);
h1 = figure;
plot(V,g_ss_h,'k--','LineWidth',4);
hold on;
plot(V,g_ss_fit,'k','LineWidth',4);
% legend('Luo and Rudy','Fitted');
xlabel('Voltage (mV)');
ylabel('h_{ss}');
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
plot(V,tau_h,'k--','LineWidth',4);
hold on;
plot(V,tau_fit,'k','LineWidth',4);
% legend('Luo and Rudy','Fitted');
xlabel('Voltage (mV)');
ylabel('\tau_h (ms)');
set(gca,'FontSize',28);
xlim([-120 60]);
set(gca,'XTick',-120:30:60);
% set(gca,'YTick',0:0.2:1);
xticklabels({-120,'',-60,'',0,'',60});
% yticklabels({0,'','','','',1});
set(gca,'LineWidth',3);
set(gca,'xgrid','on');

print_figure(h1,output_dir,'g_ss_h');
print_figure(h2,output_dir,'tau_h');

