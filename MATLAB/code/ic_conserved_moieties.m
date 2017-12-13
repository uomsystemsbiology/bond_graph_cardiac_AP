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

%% Calculate value of conserved moieties
W_i = 38;
W_e = 5.182;

perm = [4 3 39 40 41 43 45 42 44 46 37 38 6 5 59 61 63 65 60 62 64 66 67 ...
    69 71 73 68 70 72 74 2 1 47 48 49 50 51 52 53 54 55 56 57 58 25 26 ...
    27 28 29 30 21 35 36 31 32 33 34 23 24 20 19 22 18 12 13 14 15 16 17 ...
    10 11 8 9 7];

[~, X0] = cardiac_AP_exported(0,'init');
X0_IC1 = X0(perm)';

[~, X0] = cardiac_AP_exported(0,'init');
X0(4) = X0(4) - 1*W_i; % Ki
X0(6) = X0(6) + 1*W_i; % Nai
X0_IC2 = X0(perm)';

[~, X0] = cardiac_AP_exported(0,'init');
X0(3) = X0(3) + 1*W_e; % x_Ke
X0(4) = X0(4) - 1*W_e; % x_Ki
X0(5) = X0(5) - 1*W_e; % x_Nae
X0(6) = X0(6) + 1*W_e; % x_Nai
X0_IC3 = X0(perm)';

idx_stat_e = [2 14 32 60:63]; % Extracellular ions constant
idx_stat_ions = [1:2 13:14 31:32 60:63]; % Ions constant
idx_stat_var = 60:63; % All ions variable

load([storage_dir 'cardiac_AP_cm.mat']);

F = 96485;
G_var_ions(:,end) = G_var_ions(:,end)/F;
G_var_int(:,end) = G_var_int(:,end)/F;

cm_var_1 = G_var_ions*X0_IC1(setdiff(1:74,idx_stat_var));
cm_var_2 = G_var_ions*X0_IC2(setdiff(1:74,idx_stat_var));
cm_var_3 = G_var_ions*X0_IC3(setdiff(1:74,idx_stat_var));

cm_e_1 = G_var_int*X0_IC1(setdiff(1:74,idx_stat_e));
cm_e_2 = G_var_int*X0_IC2(setdiff(1:74,idx_stat_e));
cm_e_3 = G_var_int*X0_IC3(setdiff(1:74,idx_stat_e));