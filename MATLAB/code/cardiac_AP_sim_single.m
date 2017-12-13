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

%% Options
save_figures = true;

%% Define constants
R = 8.314;
T = 310;
F = 96485;

N_A = 6.022e23;
C_m = 153400; % unit fF

%% Define volumes (unit pL)
W_i = 38;
W_e = 5.182;

%% Run simulation
tspan = [0 0.8];
[VOI,STATES,CONSTANTS] = cardiac_AP_exported(tspan,'con');
% Compute algebraic variables
[RATES,ALGEBRAIC] = cardiac_AP_exported_computeRates(VOI, STATES, CONSTANTS);
ALGEBRAIC = cardiac_AP_exported_computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);
t = VOI;

%% Plot membrane voltage
h1 = figure;
V = 1000*ALGEBRAIC(:,15);
plot(1000*(t-0.3),V,'k','LineWidth',2);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
xlim(1000*[-0.1 0.4]);
set(gca,'FontSize',16);

%% Plot intracellular calcium concentration
Cai = STATES(:,2)/W_i;

h2 = figure;
plot(1000*(t-0.3),1e3*Cai,'k','LineWidth',2);
xlabel('Time (sec)','FontSize',16);
ylabel('[Ca^{2+}]_i ({\mu}M)','FontSize',16);
xlim(1000*[-0.1 0.4]);
set(gca,'FontSize',16);

%% Plot ion channel currents
I_Na = F*ALGEBRAIC(:,938)/C_m;
I_K1 = F*ALGEBRAIC(:,1014)/C_m;
I_K = F*ALGEBRAIC(:,967)/C_m;
I_Kp = F*ALGEBRAIC(:,989)/C_m;
I_Ca = 2*F*(ALGEBRAIC(:,1038)+ALGEBRAIC(:,1057))/C_m;
I_K_LCC = F*(ALGEBRAIC(:,1082)+ALGEBRAIC(:,1102))/C_m;

h3 = figure;
hold on;
plot(1000*(t-0.3),I_Na,'LineWidth',2);
plot(1000*(t-0.3),I_K1,'LineWidth',2);
plot(1000*(t-0.3),I_K,'LineWidth',2);
plot(1000*(t-0.3),I_Kp,'LineWidth',2);
plot(1000*(t-0.3),I_Ca,'LineWidth',2);
plot(1000*(t-0.3),I_K_LCC,'LineWidth',2);

xlabel('Time (ms)');
ylabel('Current ({\mu}A/{\mu}F)');
legend('I_{Na}','I_{K1}','I_{K}','I_{Kp}','I_{Ca,L}','I_{K,L}');
xlim(1000*[-0.1 0.4]);
ylim([-10 6]);
set(gca,'FontSize',16);
box on;

%% Plot transporter and gating currents
I_NaK = -ALGEBRAIC(:,931)/C_m; % unit uA/uF
I_NCX = -ALGEBRAIC(:,770)/C_m; % unit uA/uF

I_Na_gate = F*ALGEBRAIC(:,2180);
I_K1_gate = F*ALGEBRAIC(:,1166);
I_K_gate = F*ALGEBRAIC(:,1507);
I_Kp_gate = F*ALGEBRAIC(:,1144);
I_LCC_gate = F*ALGEBRAIC(:,1823);
I_gate = (I_Na_gate + I_K1_gate + I_K_gate + I_Kp_gate + I_LCC_gate)/C_m;

h4 = figure;
hold on;
plot(1000*(t-0.3),I_NaK,'LineWidth',2);
plot(1000*(t-0.3),I_NCX,'LineWidth',2);
plot(1000*(t-0.3),I_gate,'LineWidth',2);

xlabel('Time (ms)');
ylabel('Current ({\mu}A/{\mu}F)');
legend('I_{NaK}','I_{NCX}','I_{gate}');
xlim(1000*[-0.1 0.4]);
ylim([-1.5 1.5]);
set(gca,'FontSize',16);
box on;


%% Plot power dissipation
[LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = cardiac_AP_exported_createLegends();

Idx_Af = find(LEGEND_ALGEBRAIC(:,1) == 'A' & LEGEND_ALGEBRAIC(:,2) == 'f' & ...
    LEGEND_ALGEBRAIC(:,3) == '_');

array_Re = [];
for i_Re = 1:length(Idx_Af)
    idx_Re = Idx_Af(i_Re);
    str_Af_full = LEGEND_ALGEBRAIC(idx_Re,:);
    idx_space = find(str_Af_full == ' ',1);
    str_Re_name = str_Af_full(4:idx_space-1);
    length_Re_name = length(str_Re_name);

    array_Re = [array_Re; idx_Re ...
        find(all((LEGEND_ALGEBRAIC(:,1:4+length_Re_name) == ['Ar_' str_Re_name ' '])')) ...
        find(all((LEGEND_ALGEBRAIC(:,1:3+length_Re_name) == ['v_' str_Re_name ' '])'))];
end

P = zeros(length(t),1);
for i_Re = 1:size(array_Re,1)
    P = P + (ALGEBRAIC(:,array_Re(i_Re,1))-ALGEBRAIC(:,array_Re(i_Re,2))) .* ALGEBRAIC(:,array_Re(i_Re,3));
end

h5 = figure;
plot(1000*(t-0.3),P/1e3,'k','LineWidth',2);
xlabel('Time (ms)');
ylabel('Power (pW)');
xlim(1000*[-0.1 0.4]);
set(gca,'FontSize',16);

Idx_post_stim = t>0.3;
num_pre_stim = sum(t<=0.3);
E = cumtrapz(t(Idx_post_stim),P(Idx_post_stim));

t_beat = t-0.3;
V_peak = max(V);
V_end = V(end);
V_90 = 0.9*V_end + 0.1*V_peak;
Idx_repolarised = find((V <= V_90) & (t_beat > 0.02),1);
Idx_cross = [Idx_repolarised-1 Idx_repolarised];
APD = interp1(V(Idx_cross),t_beat(Idx_cross),V_90);

h6 = figure;
plot(1000*[-0.1; 0; t(Idx_post_stim)-0.3],[0; 0; E/1e3],'k','LineWidth',2);
hold on;
plot(1000*t_beat(Idx_repolarised),E(Idx_repolarised-num_pre_stim)/1e3,'k.','MarkerSize',30);
text(1000*t_beat(Idx_repolarised)-10,E(Idx_repolarised-num_pre_stim)/1e3-3,['E = ' num2str(round(E(Idx_repolarised-num_pre_stim)/1e3,1)) ' pJ'],...
    'Color','k','FontSize',16);
xlabel('Time (ms)');
ylabel('Energy (pJ)');
xlim(1000*[-0.1 0.4]);
set(gca,'FontSize',16);

%% Print figures
if save_figures
    print_figure(h1,output_dir,'cardiac_AP_Vm');
    print_figure(h2,output_dir,'cardiac_AP_Cai');
    print_figure(h3,output_dir,'cardiac_AP_ion_currents');
    print_figure(h4,output_dir,'cardiac_AP_other_currents');
    print_figure(h5,output_dir,'cardiac_AP_power');
    print_figure(h6,output_dir,'cardiac_AP_energy');
end
