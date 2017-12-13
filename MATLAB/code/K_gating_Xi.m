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

%% Calculate bond graph constants

A = exp(-56.26/32.1);
K = 1/(1+A*exp(-120/32.1));

alpha0 = K;
beta0 = A*K;

zf = 0;
zr = R*T/32.1/F;

V = -120:1:60;

alpha = alpha0*exp(zf*F*V/R/T);
beta = beta0*exp(zr*F*V/R/T);

g_ss = alpha./(alpha + beta);
tau = 1./(alpha + beta);

g_ss_LR = 1./(1+exp((V-56.26)/32.1));

h1 = figure;
hold on;
plot(V,g_ss,'k--','LineWidth',4);
plot(V,g_ss_LR,'k','LineWidth',4)
% legend('Luo and Rudy','Fitted');
xlabel('Voltage (mV)');
ylabel('Xi_{ss}');
set(gca,'FontSize',28);
box on;
xlim([-120 60]);
ylim([0 1]);
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
ylabel('\tau_{Xi} (ms)');
set(gca,'FontSize',28);
xlim([-120 60]);
set(gca,'XTick',-120:30:60);
% set(gca,'YTick',0:0.2:1);
xticklabels({-120,'',-60,'',0,'',60});
% yticklabels({0,'','','','',1});
set(gca,'LineWidth',3);
set(gca,'xgrid','on');

print_figure(h1,output_dir,'g_ss_Xi');
print_figure(h2,output_dir,'tau_Xi');