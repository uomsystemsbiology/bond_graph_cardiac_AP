clear;
clc;
close all;

%% Script options
run_optimisation = false;
marker_size = 30;
num_runs = 50;
swarm_size = 1000;
save_figures = true;

%% Set up directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];
dataset_func_dir = [code_dir 'dataset_functions' filesep];

addpath(dataset_func_dir);

%% Load Kimura et al. extracellular sodium data
file_path = [data_dir 'Kimura_Nae.csv'];
array_data = csvread(file_path,1,0);

V = transpose((-110:10:50)); % unit V
num_V = length(V);
Nae = [17.5; 35; 70; 140];
num_Nae = length(Nae);
array_I = zeros(num_Nae,num_V);

V_data = array_data(:,1);
[~,idx_unique] = unique(V_data);
V_data = V_data(idx_unique);

for i_Nae = 1:num_Nae
    I_data = array_data(idx_unique,1+i_Nae);
    array_I(i_Nae,:) = interp1(V_data,I_data,V);
end

struct_Kimura_Nae = struct('V',V/1000,...
    'Nae',Nae,...
    'array_I',array_I);

%% Load Kimura et al. extracellular calcium data
file_path = [data_dir 'Kimura_Cae.csv'];
array_data = csvread(file_path,1,0);

V = transpose((-110:10:40)); % unit V
num_V = length(V);
Cae = [1; 4];
num_Cae = length(Cae);
array_I = zeros(num_Cae,num_V);

V_data = array_data(:,1);
[~,idx_unique] = unique(V_data);
V_data = V_data(idx_unique);

for i_Cae = 1:num_Cae
    I_data = array_data(idx_unique,1+i_Cae);
    array_I(i_Cae,:) = interp1(V_data,I_data,V);
end

struct_Kimura_Cae = struct('V',V/1000,...
    'Cae',Cae,...
    'array_I',array_I);

%% Load Beuckelmann and Wier data
file_path = [data_dir 'BW_Cai.csv'];
array_data = csvread(file_path,1,0);

V = array_data(:,1); % unit V
array_I = array_data(:,2);

struct_BW = struct('V',V/1000,...
    'array_I',array_I);

%% Compile all data sources into a struct
struct_NCX_data = struct('Kimura_Nae',struct_Kimura_Nae,...
    'Kimura_Cae',struct_Kimura_Cae,...
    'BW',struct_BW);

%% Formulate minimisation problem for Nakao and Gadsby data
error_func = @(params_vec)NCX_total_error(params_vec,struct_NCX_data);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [zeros(10,1); -Inf];
ub = inf(11,1);

options_unc = optimoptions('fminunc','MaxFunEvals',10000);
options_ps = optimoptions('particleswarm','UseParallel',true,'SwarmSize',swarm_size);

%% Run parameter optimisation
if run_optimisation
    fmin = 1e10;
    tic
    for i_runs = 1:num_runs
        i_runs
        [params_vec,fval,exitflag,output] = particleswarm(error_func,11,lb,ub,options_ps);
        [params_vec,fval,exitflag,output,grad,hessian] = fminunc(error_func,params_vec,options_unc);
        toc
        if fval < fmin
            params_vec_min = params_vec;
            fmin = fval
            output_min = output;
            hessian_min = hessian;
        end
    end
    params_vec = params_vec_min;
    fval = fmin;
    output = output_min;
    hessian = hessian_min;
    
    % Scale kappa to achieve realistic cycling velocity
    struct_input = struct('Nai',0,...
        'Cae',1,...
        'Cai',430e-6,...
        'V',-0.11,...
        'Nae',140);
    scaling_ratio = 700/NCX_vss_fitting(params_vec,struct_input,1);
    params_vec(1:2) = scaling_ratio*params_vec(1:2);
    % Save parameters
    save([storage_dir 'NCX_fitting_results.mat'],'params_vec','options_ps','options_unc','hessian');
else
    load([storage_dir 'NCX_fitting_results.mat']);
end

%% Plot comparison to Kimura Nae data
% Experimental conditions:
struct_input = struct('Nai',0,...
    'Cae',1,...
    'Cai',430e-6,...
    'V',struct_Kimura_Nae.V(1),...
    'Nae',struct_Kimura_Nae.Nae(4));

Nae_vec = struct_Kimura_Nae.Nae;
V_vec = struct_Kimura_Nae.V;

% Normalising factors
normalising_factor_data = abs(struct_Kimura_Nae.array_I(4,1));
normalising_factor_model = -abs(NCX_vss_fitting(params_vec,struct_input,1));

r_I_model = zeros(length(Nae_vec),length(V_vec));
r_I_data = zeros(length(Nae_vec),length(V_vec));
for i_Nae = 1:length(Nae_vec)
    struct_input.Nae = Nae_vec(i_Nae);
    for i_V = 1:length(V_vec)
        struct_input.V = V_vec(i_V);
        r_I_model(i_Nae,i_V) = NCX_vss_fitting(params_vec,struct_input,1)/normalising_factor_model;
        r_I_data(i_Nae,i_V) = struct_Kimura_Nae.array_I(i_Nae,i_V)/normalising_factor_data;
    end
end

h1 = figure;
hold on;
box on;
plot(1000*V_vec,r_I_data,'.','MarkerSize',marker_size);
ax = gca;
ax.ColorOrderIndex = 1;
h_line = plot(1000*V_vec,r_I_model,'LineWidth',2);
xlim([-120 60]);
ylim([-1.2 0.2]);
xlabel('Voltage (mV)');
ylabel('Normalised cycling velocity');
set(gca,'FontSize',16);
h_leg = legend(h_line,{'[Na^+]_e = 17.5 mM', '[Na^+]_e = 35 mM', '[Na^+]_e = 70 mM', '[Na^+]_e = 150 mM'});
h_leg.Location = 'southeast';
h_leg.FontSize = 14;

%% Plot comparison to Kimura Cae data
% Experimental conditions:
struct_input = struct('Nai',10,...
    'Cae',struct_Kimura_Cae.Cae(end),...
    'Cai',172e-6,...
    'V',struct_Kimura_Cae.V(end),...
    'Nae',140);

Cae_vec = struct_Kimura_Cae.Cae;
V_vec = struct_Kimura_Cae.V;

% Normalising factors
normalising_factor_data = abs(struct_Kimura_Cae.array_I(end,end));
normalising_factor_model = -abs(NCX_vss_fitting(params_vec,struct_input,1));

r_I_model = zeros(length(Cae_vec),length(V_vec));
r_I_data = zeros(length(Cae_vec),length(V_vec));
for i_Cae = 1:length(struct_Kimura_Cae.Cae)
    struct_input.Cae = struct_Kimura_Cae.Cae(i_Cae);
    for i_V = 1:length(struct_Kimura_Cae.V)
        struct_input.V = struct_Kimura_Cae.V(i_V);
        r_I_model(i_Cae,i_V) = NCX_vss_fitting(params_vec,struct_input,1)/normalising_factor_model;
        r_I_data(i_Cae,i_V) = struct_Kimura_Cae.array_I(i_Cae,i_V)/normalising_factor_data;
    end
end

h2 = figure;
hold on;
box on;
plot(1000*V_vec,r_I_data,'.','MarkerSize',marker_size);
ax = gca;
ax.ColorOrderIndex = 1;
h_line = plot(1000*V_vec,r_I_model,'LineWidth',2);
xlim([-120 60]);
ylim([-0.6 1.2]);
xlabel('Voltage (mV)');
ylabel('Normalised cycling velocity');
set(gca,'FontSize',16);
h_leg = legend(h_line,{'[Ca^{2+}]_e = 1 mM','[Ca^{2+}]_e = 4 mM'});
h_leg.Location = 'northwest';
h_leg.FontSize = 14;

%% Plot comparison to Beuckelmann and Wier data
% Experimental conditions:
struct_input = struct('Nai',15,...
    'Cae',2,...
    'Cai',450e-6,...
    'V',struct_BW.V(end),...
    'Nae',135);

V_vec = struct_BW.V;

% Normalising factors
normalising_factor_data = abs(struct_BW.array_I(end));
normalising_factor_model = -abs(NCX_vss_fitting(params_vec,struct_input,1));

r_I_model = zeros(1,length(V_vec));
r_I_data = zeros(1,length(V_vec));
for i_V = 1:length(struct_BW.V)
    struct_input.V = struct_BW.V(i_V);
    r_I_model(i_V) = NCX_vss_fitting(params_vec,struct_input,1)/normalising_factor_model;
    r_I_data(i_V) = struct_BW.array_I(i_V)/normalising_factor_data;
end

h3 = figure;
hold on;
box on;
plot(1000*V_vec,r_I_data,'k.','MarkerSize',marker_size);
plot(1000*V_vec,r_I_model,'k','LineWidth',2);
xlim([-180 80]);
ylim([-1.1 1.1]);
xlabel('Voltage (mV)');
ylabel('Normalised cycling velocity');
set(gca,'FontSize',16);

%% Print figures
if save_figures
    print_figure(h1,output_dir,'NCX_Nae_fit');
    print_figure(h2,output_dir,'NCX_Cae_fit');
    print_figure(h3,output_dir,'NCX_Cai_fit');
end

