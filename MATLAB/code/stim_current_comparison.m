clear;
clc;
close all;

%% Set up directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];
save_figures = true;

%% Define constants
R = 8.314; % unit J/mol/K
T = 310;
F = 96485;

font_size = 14;

%% Indices for conserved moieties
perm = [4 3 39 40 41 43 45 42 44 46 37 38 6 5 59 61 63 65 60 62 64 66 67 ...
    69 71 73 68 70 72 74 2 1 47 48 49 50 51 52 53 54 55 56 57 58 25 26 ...
    27 28 29 30 21 35 36 31 32 33 34 23 24 20 19 22 18 12 13 14 15 16 17 ...
    10 11 8 9 7];
iperm(perm) = 1:length(perm);
idx_stat_e = [2 14 32 60:63]; % Extracellular ions constant
idx_stat_ions = [1:2 13:14 31:32 60:63]; % Ions constant
idx_stat_var = 60:63; % All ions variable

load([storage_dir 'cardiac_AP_cm.mat']);

G_charge = G_var_int(end,:);
G_charge(end) = G_charge(end)/F;
for i_stat = idx_stat_e
    G_charge = [G_charge(1:i_stat-1) 0 G_charge(i_stat:end)];
end
G_charge = G_charge(iperm);

%% Plot effects of stimulus current on long-term behaviour
struct_var_con = extract_results([storage_dir 'conservative_stim_result_var' filesep 'AP_sim_states_con_'],30,G_charge);
struct_var_ncs = extract_results([storage_dir 'nonconservative_stim_result_var' filesep 'AP_sim_states_con_'],30,G_charge);

struct_var_con.APD = calc_APD_array(struct_var_con.t,struct_var_con.V);
struct_var_ncs.APD = calc_APD_array(struct_var_ncs.t,struct_var_ncs.V);

struct_var_con.peak_Ca = calc_peak_Ca(struct_var_con.t,struct_var_con.Cai);
struct_var_ncs.peak_Ca = calc_peak_Ca(struct_var_ncs.t,struct_var_ncs.Cai);

h_fig = figure;
subplot(6,3,1);
plot(struct_var_con.t(301:1000:end)/60,1000*struct_var_con.V(301:1000:end),'b','LineWidth',2);
hold on;
plot(struct_var_ncs.t(301:1000:end)/60,1000*struct_var_ncs.V(301:1000:end),'r','LineWidth',2);
xlim([0 30]);
ylabel('V_{dia} (mV)');
set(gca,'FontSize',font_size);
title('A','FontSize',28);

subplot(6,3,4);
plot(struct_var_con.t(1:1000:end-1)/60,struct_var_con.APD,'b','LineWidth',2);
hold on;
plot(struct_var_ncs.t(1:1000:end-1)/60,struct_var_ncs.APD,'r','LineWidth',2);
xlim([0 30]);
ylabel('APD (ms)');
set(gca,'FontSize',font_size);

subplot(6,3,7);
plot(struct_var_con.t/60,struct_var_con.Ki,'b','LineWidth',2);
hold on;
plot(struct_var_ncs.t/60,struct_var_ncs.Ki,'r','LineWidth',2);
xlim([0 30]);
ylabel('[K^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(6,3,10);
plot(struct_var_con.t/60,struct_var_con.Nai,'b','LineWidth',2);
hold on;
plot(struct_var_ncs.t/60,struct_var_ncs.Nai,'r','LineWidth',2);
xlim([0 30]);
ylabel('[Na^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(6,3,13);
plot(struct_var_con.t(1:1000:end-1)/60,1000*struct_var_con.peak_Ca,'b','LineWidth',2);
hold on;
plot(struct_var_ncs.t(1:1000:end-1)/60,1000*struct_var_ncs.peak_Ca,'r','LineWidth',2);
xlim([0 30]);
ylabel('Peak [Ca^{2+}]_i ({\mu}M)');
xlabel('Time (mins)');
set(gca,'FontSize',font_size);

initial_charge = struct_var_con.charge(1);

subplot(6,3,16);
plot(struct_var_con.t(1:3001),struct_var_con.charge(1:3001)-initial_charge,'b','LineWidth',2);
hold on;
plot(struct_var_ncs.t(1:3001),struct_var_ncs.charge(1:3001)-initial_charge,'r','LineWidth',2);
ylabel('{\Delta}Charge (fmol)');
xlabel('Time (sec)');
set(gca,'FontSize',font_size);

clear struct_var_con struct_var_ncs;

struct_e_con = extract_results([storage_dir 'conservative_stim_result' filesep 'AP_sim_states_con_'],30,G_charge);
struct_e_ncs = extract_results([storage_dir 'nonconservative_stim_result' filesep 'AP_sim_states_ncs_'],30,G_charge);

struct_e_con.APD = calc_APD_array(struct_e_con.t,struct_e_con.V);
struct_e_ncs.APD = calc_APD_array(struct_e_ncs.t,struct_e_ncs.V);

struct_e_con.peak_Ca = calc_peak_Ca(struct_e_con.t,struct_e_con.Cai);
struct_e_ncs.peak_Ca = calc_peak_Ca(struct_e_ncs.t,struct_e_ncs.Cai);

subplot(6,3,2);
plot(struct_e_con.t(301:1000:end)/60,1000*struct_e_con.V(301:1000:end),'b','LineWidth',2);
hold on;
plot(struct_e_ncs.t(301:1000:end)/60,1000*struct_e_ncs.V(301:1000:end),'r','LineWidth',2);
xlim([0 30]);
ylabel('V_{dia} (mV)');
set(gca,'FontSize',font_size);
title('B','FontSize',28);

subplot(6,3,5);
plot(struct_e_con.t(1:1000:end-1)/60,struct_e_con.APD,'b','LineWidth',2);
hold on;
plot(struct_e_ncs.t(1:1000:end-1)/60,struct_e_ncs.APD,'r','LineWidth',2);
xlim([0 30]);
ylabel('APD (ms)');
set(gca,'FontSize',font_size);

subplot(6,3,8);
plot(struct_e_con.t/60,struct_e_con.Ki,'b','LineWidth',2);
hold on;
plot(struct_e_ncs.t/60,struct_e_ncs.Ki,'r','LineWidth',2);
xlim([0 30]);
ylabel('[K^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(6,3,11);
plot(struct_e_con.t/60,struct_e_con.Nai,'b','LineWidth',2);
hold on;
plot(struct_e_ncs.t/60,struct_e_ncs.Nai,'r','LineWidth',2);
xlim([0 30]);
ylabel('[Na^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(6,3,14);
plot(struct_e_con.t(1:1000:end-1)/60,1000*struct_e_con.peak_Ca,'b','LineWidth',2);
hold on;
plot(struct_e_ncs.t(1:1000:end-1)/60,1000*struct_e_ncs.peak_Ca,'r','LineWidth',2);
xlim([0 30]);
ylabel('Peak [Ca^{2+}]_i ({\mu}M)');
xlabel('Time (mins)');
set(gca,'FontSize',font_size);

subplot(6,3,17);
plot(struct_e_con.t(1:3001),struct_e_con.charge(1:3001)-initial_charge,'b','LineWidth',2);
hold on;
plot(struct_e_ncs.t(1:3001),struct_e_ncs.charge(1:3001)-initial_charge,'r','LineWidth',2);
ylabel('{\Delta}Charge (fmol)');
xlabel('Time (sec)');
set(gca,'FontSize',font_size);

clear struct_e_con struct_e_ncs;

struct_ie_con = extract_results([storage_dir 'ie_stat_conservative_stim_result' filesep 'AP_sim_states_'],30,G_charge);
struct_ie_ncs = extract_results([storage_dir 'ie_stat_nonconservative_stim_result' filesep 'AP_sim_states_'],30,G_charge);

struct_ie_con.APD = calc_APD_array(struct_ie_con.t,struct_ie_con.V);
struct_ie_ncs.APD = calc_APD_array(struct_ie_ncs.t,struct_ie_ncs.V);

subplot(6,3,3);
plot(struct_ie_con.t(301:1000:end)/60,1000*struct_ie_con.V(301:1000:end),'b','LineWidth',2);
hold on;
plot(struct_ie_ncs.t(301:1000:end)/60,1000*struct_ie_ncs.V(301:1000:end),'r','LineWidth',2);
xlim([0 30]);
ylabel('V_{dia} (mV)');
title('C');
set(gca,'FontSize',font_size);
ylim([-89 -87]);
legend('Conservative','Nonconservative');
title('C','FontSize',28);

subplot(6,3,6);
plot(struct_ie_con.t(1:1000:end-1)/60,struct_ie_con.APD,'b','LineWidth',2);
hold on;
plot(struct_ie_ncs.t(1:1000:end-1)/60,struct_ie_ncs.APD,'r','LineWidth',2);
xlim([0 30]);
ylabel('APD (ms)');
set(gca,'FontSize',font_size);

subplot(6,3,9);
plot(struct_ie_con.t/60,struct_ie_con.Ki,'b','LineWidth',2);
hold on;
plot(struct_ie_ncs.t/60,struct_ie_ncs.Ki,'r','LineWidth',2);
xlim([0 30]);
ylabel('[K^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(6,3,12);
plot(struct_ie_con.t/60,struct_ie_con.Nai,'b','LineWidth',2);
hold on;
plot(struct_ie_ncs.t/60,struct_ie_ncs.Nai,'r','LineWidth',2);
xlim([0 30]);
ylabel('[Na^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(6,3,15);
plot(struct_ie_con.t/60,1000*struct_ie_con.Cai,'b','LineWidth',2);
hold on;
plot(struct_ie_ncs.t/60,1000*struct_ie_ncs.Cai,'r','LineWidth',2);
xlim([0 30]);
ylabel('[Ca^{2+}]_i ({\mu}M)');
xlabel('Time (mins)');
set(gca,'FontSize',font_size);
ylim([0 1]);

subplot(6,3,18);
plot(struct_ie_con.t(1:3001),struct_ie_con.charge(1:3001)-initial_charge,'b','LineWidth',2);
hold on;
plot(struct_ie_ncs.t(1:3001),struct_ie_ncs.charge(1:3001)-initial_charge,'r','LineWidth',2);
ylabel('{\Delta}Charge (fmol)');
xlabel('Time (sec)');
set(gca,'FontSize',font_size);

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 30 40]);
print_figure(h_fig,output_dir,'stim_current_comparison');
