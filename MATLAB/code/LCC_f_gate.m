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
run_optimisation = false;

%% Find parameters for f-gate
V = transpose(-90:1:50);

f_ss = 1./(1+exp((V+32)/8)) + 0.6./(1+exp((50-V)/20)); 
tau_f = 1./(0.0197*exp(-(0.0337*(V+10)).^2)+0.02); % Unit ms

error_func = @(params) f_gate_error(params,V,f_ss,tau_f,false);

lb = [1e-6; -5; 1e-6; -0.5; 1e-6; -10; 1e-6; -10];
ub = [1; 10; 1; 10; 1; 5; 1; 1];
options_unc = optimoptions('fminunc','MaxFunEvals',10000);
options_ps = optimoptions('particleswarm','HybridFcn',@fmincon,...
    'Display','iter',...
    'FunctionTolerance',5e-3);

if run_optimisation
    [params_vec,fval,exitflag,output] = particleswarm(error_func,8,lb,ub,options_ps);
	save([storage_dir 'LCC_f_parameters.mat'],'params_vec');
else
    load([storage_dir 'LCC_f_parameters.mat']);
end

V_ext = transpose(-120:1:60);
f_ss_ext = 1./(1+exp((V_ext+32)/8)) + 0.6./(1+exp((50-V_ext)/20));
tau_f_ext = 1./(0.0197*exp(-(0.0337*(V_ext+10)).^2)+0.02); % Unit ms

% Plot a comparison between Luo-Rudy and BG f-gates in simulations run in
% fitting procedure
%f_gate_error(params_vec,V_ext,f_ss_ext,tau_f_ext,true);

[f_ss_fit,tau_fit] = f_gate_func(params_vec,V_ext/1000);

h1 = figure;
hold on;
plot(V_ext,f_ss_ext,'k:','LineWidth',3);
plot(V_ext,f_ss_fit,'k','LineWidth',3);
xlabel('Voltage (mV)');
ylabel('f_{ss}');
xlim([-120 60]);
ylim([0 1]);
legend('LRd','BG')
set(gca,'FontSize',24);
box on;
set(gca,'XTick',-120:30:60);
xticklabels({-120,'',-60,'',0,'',60});

print_figure(h1,output_dir,'g_ss_f');

%% Plot comparison of gates with an AP-like voltage input
CellML_data = csvread([data_dir 'luo_rudy_1994_data.csv'],1,0);

t_LRd = CellML_data(:,1); % unit ms
V_LRd = CellML_data(:,5); % unit mV

f_0 = 1;
tspan = [0 t_LRd(end)];

dfdt_LRd = @(t,f) df_LRd(f,t,t_LRd,V_LRd);

options = odeset('RelTol',1e-3,'AbsTol',0.0005);
[t1,f] = ode15s(dfdt_LRd,tspan,f_0,options);

dxdt_BG = @(t,x) df_BG(x,t,t_LRd,V_LRd/1000,params_vec);
[t2,x] = ode15s(dxdt_BG,tspan,[1;0;0],options);

h2 = figure;
plot(t1,f,'k:','LineWidth',3);
hold on;
plot(t2,x(:,1)+x(:,3),'k','LineWidth',3);
xlabel('Time (ms)');
ylabel('f');
legend('LRd','BG','Location','southeast')
set(gca,'FontSize',24)

h3 = figure;
plot(t_LRd,V_LRd,'k','LineWidth',3);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
set(gca,'FontSize',24)

print_figure(h2,output_dir,'f_AP_input');
print_figure(h3,output_dir,'AP_input');

