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

%% Plot I-V curves
cNao_st = 140;
cNai_st = 10;

A_geo = 1.534e-4; % Unit cm^2
G_Na = 16*A_geo; % Unit mA/V

V = (-120:1:60)/1000;
E_Na = R*T/F*log(cNao_st/cNai_st);
I_lin = G_Na*(V-E_Na);

E_Na_norm = F*E_Na/(R*T);
G_GHK = 2*G_Na*(1-exp(E_Na_norm))/(cNai_st - cNao_st*exp(E_Na_norm))*R*T/F; % Unit mA/mM
V_norm = F*V/R/T;
GHK_factor = V_norm./(1-exp(-V_norm));
idx_nan = isnan(GHK_factor);
GHK_factor(idx_nan) = 1;
I_GHK = G_GHK*GHK_factor.*(cNai_st - cNao_st*exp(-V_norm));
% I_GHK = G_GHK*GHK_factor.*(100 - cNao_st*exp(-V_norm));

h = figure;
plot(1000*V,1e6*I_lin,'k--',1000*V,1e6*I_GHK,'k','LineWidth',2);
legend('LRd','BG','Location','southeast');
ylabel('Current (nA)');
xlabel('Voltage (mV)');
set(gca,'FontSize',16);

print_figure(h,output_dir,'Na_IV_curve');