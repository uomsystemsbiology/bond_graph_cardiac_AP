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

%% Define volumes (unit pL)
W_i = 38;
W_e = 5.182;

%% Calculate gate transition parameters
% Calculate parameters for X_i gate
A = exp(7.488/5.98);
K = 1/(1+A*exp(-60/5.98));

alpha0 = K; % Unit ms^-1
beta0 = A*K; % Unit ms^-1

zf = 0;
zr = -R*T/5.98/F;

% Plot fit of parameters to steady-state curve
V = -120:1:60;

alpha = alpha0*exp(zf*F*V/R/T);
beta = beta0*exp(zr*F*V/R/T);

g_ss = alpha./(alpha + beta);
tau = 1./(alpha + beta);

g_ss_LR = 1./(1+exp((7.488-V)/5.98));

h1 = figure;
hold on;
plot(V,g_ss,'k','LineWidth',4);
plot(V,g_ss_LR,'k--','LineWidth',4)
% legend('Luo and Rudy','Fitted','Location','northwest');
xlabel('Voltage (mV)');
ylabel('Kp_{ss}');
set(gca,'FontSize',28);
box on;
xlim([-120 60]);
set(gca,'XTick',-120:30:60);
set(gca,'YTick',0:0.2:1);
xticklabels({-120,'',-60,'',0,'',60});
yticklabels({0,'','','','',1});
set(gca,'LineWidth',3);
grid on;

colorOrder = get(gca, 'ColorOrder');

h2 = figure;
plot(V,tau,'k','LineWidth',4);
xlabel('Voltage (mV)');
ylabel('\tau_{Kp} (ms)');
set(gca,'FontSize',28);
xlim([-120 60]);
set(gca,'XTick',-120:30:60);
% set(gca,'YTick',0:0.2:1);
xticklabels({-120,'',-60,'',0,'',60});
% yticklabels({0,'','','','',1});
set(gca,'LineWidth',3);
set(gca,'xgrid','on');

params_Kp = [alpha0; zf*1e3; beta0; zr*1e3];
