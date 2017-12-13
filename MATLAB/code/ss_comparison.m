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

%% Plot effects of stimulus current on long-term behaviour
struct_var_ic1 = extract_results([storage_dir 'conservative_stim_result_var' filesep 'AP_sim_states_con_'],10);
struct_var_ic2 = extract_results([storage_dir 'ic2_var' filesep 'AP_sim_states_'],10);
struct_var_ic3 = extract_results([storage_dir 'ic3_var' filesep 'AP_sim_states_'],10);

struct_var_ic1.APD = calc_APD_array(struct_var_ic1.t,struct_var_ic1.V);
struct_var_ic2.APD = calc_APD_array(struct_var_ic2.t,struct_var_ic2.V);
struct_var_ic3.APD = calc_APD_array(struct_var_ic3.t,struct_var_ic3.V);

struct_var_ic1.peak_Ca = calc_peak_Ca(struct_var_ic1.t,struct_var_ic1.Cai);
struct_var_ic2.peak_Ca = calc_peak_Ca(struct_var_ic2.t,struct_var_ic2.Cai);
struct_var_ic3.peak_Ca = calc_peak_Ca(struct_var_ic3.t,struct_var_ic3.Cai);

h_fig = figure;
subplot(5,3,1);
plot(struct_var_ic1.t(301:1000:end)/60,1000*struct_var_ic1.V(301:1000:end),'k','LineWidth',2);
hold on;
plot(struct_var_ic2.t(301:1000:end)/60,1000*struct_var_ic2.V(301:1000:end),'b','LineWidth',2);
plot(struct_var_ic3.t(301:1000:end)/60,1000*struct_var_ic3.V(301:1000:end),'r','LineWidth',2);
xlim([0 10]);
ylabel('V_{dia} (mV)');
set(gca,'FontSize',font_size);
title('A','FontSize',28);

subplot(5,3,4);
plot(struct_var_ic1.t(1:1000:end-1)/60,struct_var_ic1.APD,'k','LineWidth',2);
hold on;
plot(struct_var_ic2.t(1:1000:end-1)/60,struct_var_ic2.APD,'b','LineWidth',2);
plot(struct_var_ic3.t(1:1000:end-1)/60,struct_var_ic3.APD,'r','LineWidth',2);
xlim([0 10]);
ylabel('APD (ms)');
set(gca,'FontSize',font_size);

subplot(5,3,7);
plot(struct_var_ic1.t/60,struct_var_ic1.Ki,'k','LineWidth',2);
hold on;
plot(struct_var_ic2.t/60,struct_var_ic2.Ki,'b','LineWidth',2);
plot(struct_var_ic3.t/60,struct_var_ic3.Ki,'r','LineWidth',2);
xlim([0 10]);
ylabel('[K^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(5,3,10);
plot(struct_var_ic1.t/60,struct_var_ic1.Nai,'k','LineWidth',2);
hold on;
plot(struct_var_ic2.t/60,struct_var_ic2.Nai,'b','LineWidth',2);
plot(struct_var_ic3.t/60,struct_var_ic3.Nai,'r','LineWidth',2);
xlim([0 10]);
ylabel('[Na^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(5,3,13);
plot(struct_var_ic1.t(1:1000:end-1)/60,1000*struct_var_ic1.peak_Ca,'k','LineWidth',2);
hold on;
plot(struct_var_ic2.t(1:1000:end-1)/60,1000*struct_var_ic2.peak_Ca,'b','LineWidth',2);
plot(struct_var_ic3.t(1:1000:end-1)/60,1000*struct_var_ic3.peak_Ca,'r','LineWidth',2);
xlim([0 10]);
ylabel('Peak [Ca^{2+}]_i ({\mu}M)');
xlabel('Time (mins)');
set(gca,'FontSize',font_size);

clear struct_var_ic1 struct_var_ic2;

struct_e_ic1 = extract_results([storage_dir 'conservative_stim_result' filesep 'AP_sim_states_con_'],10);
struct_e_ic2 = extract_results([storage_dir 'ic2_e_const' filesep 'AP_sim_states_'],10);
struct_e_ic3 = extract_results([storage_dir 'ic3_e_const' filesep 'AP_sim_states_'],10);

struct_e_ic1.APD = calc_APD_array(struct_e_ic1.t,struct_e_ic1.V);
struct_e_ic2.APD = calc_APD_array(struct_e_ic2.t,struct_e_ic2.V);
struct_e_ic3.APD = calc_APD_array(struct_e_ic3.t,struct_e_ic3.V);

struct_e_ic1.peak_Ca = calc_peak_Ca(struct_e_ic1.t,struct_e_ic1.Cai);
struct_e_ic2.peak_Ca = calc_peak_Ca(struct_e_ic2.t,struct_e_ic2.Cai);
struct_e_ic3.peak_Ca = calc_peak_Ca(struct_e_ic3.t,struct_e_ic3.Cai);

subplot(5,3,2);
plot(struct_e_ic1.t(301:1000:end)/60,1000*struct_e_ic1.V(301:1000:end),'k','LineWidth',2);
hold on;
plot(struct_e_ic2.t(301:1000:end)/60,1000*struct_e_ic2.V(301:1000:end),'b','LineWidth',2);
plot(struct_e_ic3.t(301:1000:end)/60,1000*struct_e_ic3.V(301:1000:end),'r','LineWidth',2);
xlim([0 10]);
ylabel('V_{dia} (mV)');
set(gca,'FontSize',font_size);
title('B','FontSize',28);

subplot(5,3,5);
plot(struct_e_ic1.t(1:1000:end-1)/60,struct_e_ic1.APD,'k','LineWidth',2);
hold on;
plot(struct_e_ic2.t(1:1000:end-1)/60,struct_e_ic2.APD,'b','LineWidth',2);
plot(struct_e_ic3.t(1:1000:end-1)/60,struct_e_ic3.APD,'r','LineWidth',2);
xlim([0 10]);
ylabel('APD (ms)');
set(gca,'FontSize',font_size);

subplot(5,3,8);
plot(struct_e_ic1.t/60,struct_e_ic1.Ki,'k','LineWidth',2);
hold on;
plot(struct_e_ic2.t/60,struct_e_ic2.Ki,'b','LineWidth',2);
plot(struct_e_ic3.t/60,struct_e_ic3.Ki,'r','LineWidth',2);
xlim([0 10]);
ylabel('[K^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(5,3,11);
plot(struct_e_ic1.t/60,struct_e_ic1.Nai,'k','LineWidth',2);
hold on;
plot(struct_e_ic2.t/60,struct_e_ic2.Nai,'b','LineWidth',2);
plot(struct_e_ic3.t/60,struct_e_ic3.Nai,'r','LineWidth',2);
xlim([0 10]);
ylabel('[Na^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(5,3,14);
plot(struct_e_ic1.t(1:1000:end-1)/60,1000*struct_e_ic1.peak_Ca,'k','LineWidth',2);
hold on;
plot(struct_e_ic2.t(1:1000:end-1)/60,1000*struct_e_ic2.peak_Ca,'b','LineWidth',2);
plot(struct_e_ic3.t(1:1000:end-1)/60,1000*struct_e_ic3.peak_Ca,'r','LineWidth',2);
xlim([0 10]);
ylabel('Peak [Ca^{2+}]_i ({\mu}M)');
xlabel('Time (mins)');
set(gca,'FontSize',font_size);

clear struct_e_ic1 struct_e_ic2;

% struct_ie_con = extract_results([storage_dir 'conservative_stim_result_var' filesep 'AP_sim_states_con_'],30);
struct_ie_ic1 = extract_results([storage_dir 'ie_stat_conservative_stim_result' filesep 'AP_sim_states_'],10);
struct_ie_ic2 = extract_results([storage_dir 'ic2_ie_const' filesep 'AP_sim_states_'],10);
struct_ie_ic3 = extract_results([storage_dir 'ic3_ie_const' filesep 'AP_sim_states_'],10);

struct_ie_ic1.APD = calc_APD_array(struct_ie_ic1.t,struct_ie_ic1.V);
struct_ie_ic2.APD = calc_APD_array(struct_ie_ic2.t,struct_ie_ic2.V);
struct_ie_ic3.APD = calc_APD_array(struct_ie_ic3.t,struct_ie_ic3.V);

subplot(5,3,3);
plot(struct_ie_ic1.t(301:1000:end)/60,1000*struct_ie_ic1.V(301:1000:end),'k','LineWidth',2);
hold on;
plot(struct_ie_ic2.t(301:1000:end)/60,1000*struct_ie_ic2.V(301:1000:end),'b','LineWidth',2);
plot(struct_ie_ic3.t(301:1000:end)/60,1000*struct_ie_ic3.V(301:1000:end),'r','LineWidth',2);
xlim([0 10]);
ylabel('V_{dia} (mV)');
set(gca,'FontSize',font_size);
legend('IC1','IC2','IC3');
title('C','FontSize',28);

subplot(5,3,6);
plot(struct_ie_ic1.t(1:1000:end-1)/60,struct_ie_ic1.APD,'k','LineWidth',2);
hold on;
plot(struct_ie_ic2.t(1:1000:end-1)/60,struct_ie_ic2.APD,'b','LineWidth',2);
plot(struct_ie_ic3.t(1:1000:end-1)/60,struct_ie_ic3.APD,'r','LineWidth',2);
xlim([0 10]);
ylabel('APD (ms)');
set(gca,'FontSize',font_size);

subplot(5,3,9);
plot(struct_ie_ic1.t/60,struct_ie_ic1.Ki,'k','LineWidth',2);
hold on;
plot(struct_ie_ic2.t/60,struct_ie_ic2.Ki,'b','LineWidth',2);
plot(struct_ie_ic3.t/60,struct_ie_ic3.Ki,'r','LineWidth',2);
xlim([0 10]);
ylim([143.5 145.5]);
ylabel('[K^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(5,3,12);
plot(struct_ie_ic1.t/60,struct_ie_ic1.Nai,'k','LineWidth',2);
hold on;
plot(struct_ie_ic2.t/60,struct_ie_ic2.Nai,'b','LineWidth',2);
plot(struct_ie_ic3.t/60,struct_ie_ic3.Nai,'r','LineWidth',2);
xlim([0 10]);
ylim([9.5 11.5]);
ylabel('[Na^+]_i (mM)');
set(gca,'FontSize',font_size);

subplot(5,3,15);
plot(struct_ie_ic1.t/60,1000*struct_ie_ic1.Cai,'k','LineWidth',2);
hold on;
plot(struct_ie_ic2.t/60,1000*struct_ie_ic2.Cai,'b','LineWidth',2);
plot(struct_ie_ic3.t/60,1000*struct_ie_ic3.Cai,'r','LineWidth',2);
xlim([0 10]);
ylabel('[Ca^{2+}]_i ({\mu}M)');
xlabel('Time (mins)');
set(gca,'FontSize',font_size);
ylim([0 1]);

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 30 30]);
print_figure(h_fig,output_dir,'ss_comparison');
