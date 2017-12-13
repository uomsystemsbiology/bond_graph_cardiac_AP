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
R = 8.314;
T = 310;
F = 96485;

%% Steady-state gating parameters and time constants

V = transpose(-120:1:60);

g_ss_d = 1./(1+exp(-(V+10)/6.24));
tau_d = g_ss_d .* (1-exp(-(V+10)/6.24))./(0.035*(V+10));
idx_nan = find(isnan(tau_d));
tau_d(idx_nan) = g_ss_d(idx_nan)/0.035/6.24;

alpha_d = g_ss_d ./ tau_d;
beta_d = (1-g_ss_d) ./ tau_d;

%% Calculate bond graph parameters for model
K = 0.5/2.289;
alpha_0 = K*exp(10/12.48); % Unit ms^-1
z_f = 1e3/12.48*R*T/F;

beta_0 = K*exp(-10/12.48); % Unit ms^-1
z_r = -1e3/12.48*R*T/F;

params_vec = [alpha_0; z_f; beta_0; z_r];

g_ss_fit = p2gss(params_vec,V);
h1 = figure;
plot(V,g_ss_d,'k--','LineWidth',4);
hold on;
plot(V,g_ss_fit,'k','LineWidth',4);
% legend('Luo and Rudy','Fitted','Location','northwest');
xlabel('Voltage (mV)');
ylabel('d_{ss}');
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
plot(V,tau_d,'k--','LineWidth',4);
hold on;
plot(V,tau_fit,'k','LineWidth',4);
% legend('Luo and Rudy','Fitted','Location','northwest');
xlabel('Voltage (mV)');
ylabel('\tau_d (ms)');
set(gca,'FontSize',28);
xlim([-120 60]);
set(gca,'XTick',-120:30:60);
% set(gca,'YTick',0:0.2:1);
xticklabels({-120,'',-60,'',0,'',60});
% yticklabels({0,'','','','',1});
set(gca,'LineWidth',3);
set(gca,'xgrid','on');

print_figure(h1,output_dir,'g_ss_d');
print_figure(h2,output_dir,'tau_d');

save([storage_dir 'LCC_d_parameters.mat'],'params_vec');

