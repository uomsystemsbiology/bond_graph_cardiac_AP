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
alpha_j = zeros(length(V),1);
beta_j = zeros(length(V),1);

for i_V = 1:length(V)
    if V(i_V) >= -40
        alpha_j(i_V) = 0;
        beta_j(i_V) = 0.3*exp(-2.535e-7*V(i_V))/(1+exp(-0.1*(V(i_V)+32)));
    else
        alpha_j(i_V) = (-1.2714e5*exp(0.2444*V(i_V)) - 3.474e-5*exp(-0.04391*V(i_V))) * (V(i_V)+37.78) / (1 + exp(0.311*(V(i_V)+79.23)));
        beta_j(i_V) = 0.1212*exp(-0.01052*V(i_V))/(1 + exp(-0.1378*(V(i_V) + 40.14)));
    end
end

tau_j = 1./(alpha_j + beta_j);
g_ss_j = alpha_j./(alpha_j + beta_j);


%% Fit bond graph parameters to model
% params: [kf, zf, kr, zr];
error_func_alpha = @(params) square_error(alpha_j - calc_alpha(params,V/1000));
error_func_beta = @(params) square_error(beta_j - calc_beta(params,V/1000));
error_func_gss = @(params) square_error(g_ss_j - p2gss(params,V));
error_func_tau = @(params) square_error(tau_j- p2tau(params,V));

error_func = @(params) error_func_alpha(params) + error_func_beta(params) + error_func_gss(params) + error_func_tau(params);

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
    save([storage_dir 'Na_j_parameters.mat'],'params_vec');
else
    load([storage_dir 'Na_j_parameters.mat']);
end
% [params_vec,fval,exitflag,output,grad,hessian] = fminunc(error_func,params_vec,options_unc);

g_ss_fit = p2gss(params_vec,V);
h1 = figure;
plot(V,g_ss_j,'k--','LineWidth',4);
hold on;
plot(V,g_ss_fit,'k','LineWidth',4);
% legend('Luo and Rudy','Fitted');
xlabel('Voltage (mV)');
ylabel('j_{ss}');
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
plot(V,tau_j,'k--','LineWidth',4);
hold on;
plot(V,tau_fit,'k','LineWidth',4);
% legend('Luo and Rudy','Fitted');
xlabel('Voltage (mV)');
ylabel('\tau_j (ms)');
set(gca,'FontSize',28);
xlim([-120 60]);
set(gca,'XTick',-120:30:60);
% set(gca,'YTick',0:0.2:1);
xticklabels({-120,'',-60,'',0,'',60});
% yticklabels({0,'','','','',1});
set(gca,'LineWidth',3);
set(gca,'xgrid','on');

print_figure(h1,output_dir,'g_ss_j');
print_figure(h2,output_dir,'tau_j');
