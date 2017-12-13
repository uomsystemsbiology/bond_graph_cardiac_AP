clear;
clc;
close all;
alpha1_0 = 5;
beta1_0 = exp(3);

z_f1 = 10;
z_r1 = 100;

V = (-100:1:100)/1000;
alpha1 = alpha1_0*exp(z_f1*V);
beta1 = beta1_0*exp(z_r1*V);
% figure;
% plot(V,alpha1,V,beta1);
% ylim([0 50]);
% grid on;

alpha2_0 = 5;
beta2_0 = exp(3);

z_f2 = -10;
z_r2 = -100;

alpha2 = alpha2_0*exp(z_f2*V);
beta2 = beta2_0*exp(z_r2*V);
% figure;
% plot(V,alpha2,V,beta2);
% ylim([0 50]);
% grid on;

x_0 = [0.5; 0; 0.5];
options = odeset('NonNegative',1:3,'RelTol',1e-3,'AbsTol',1e-8);

sim_ode = @(t,x) dx_f(x,-100/1000);

tspan = [0; 5];
tic
[t,x] = ode15s(sim_ode,tspan,x_0,options);
toc

figure;
plot(t,x,'LineWidth',2);
legend('O_1','C','O_2');
grid on;

% Compare to analytical solution
A = [-beta1(1) alpha1(1) 0; beta1(1) -(alpha1(1)+alpha2(1)) beta2(1); 0 alpha2(1) -beta2(1)];
t_a = transpose(0:0.01:5);
x_a = zeros(length(t_a),3);
expAt = cell(length(t_a),1);
for i = 1:length(t_a)
    expAt{i} = expm(A*t_a(i));
end
tic
for i = 1:length(t_a)
    x_a(i,:) = expAt{i}*x_0;
end
toc
figure;
plot(t_a,x_a,'LineWidth',2);
legend('O_1','C','O_2');
grid on;

x_init1 = [1; 0; 0];
x_init2 = [0; 1; 0];
x_init3 = [0; 0; 1];

options = odeset('NonNegative',1:3,'RelTol',1e-3,'AbsTol',1e-8,'Refine',100);

V = (-100:5:100)/1000;
f_ss = zeros(length(V),1);
tau = zeros(length(V),1);

tic
t = transpose(0:0.01:5);
x = zeros(length(t_a),3);
for i = 1:length(V)
    sim_ode = @(t,x) dx_f(x,V(i));
    alpha1 = alpha1_0*exp(z_f1*V(i));
    beta1 = beta1_0*exp(z_r1*V(i));
    K1 = alpha1/beta1;

    alpha2 = alpha2_0*exp(z_f2*V(i));
    beta2 = beta2_0*exp(z_r2*V(i));
    K2 = alpha2/beta2;
    
    f_ss(i) = (K1+K2)/(1+K1+K2);
    
    A = [-beta1 alpha1 0; beta1 -(alpha1+alpha2) beta2; 0 alpha2 -beta2];
    expAt = cell(length(t_a),1);
    for i_t = 1:length(t_a)
        expAt{i_t} = expm(A*t_a(i_t));
    end
    
    for i_t = 1:length(t_a)
        x(i_t,:) = expAt{i_t}*x_init1;
    end
    deviation = max(abs(x(:,2)+f_ss(i)-1));
    
    idx_t = max(find( abs(x(:,2)+f_ss(i)-1) > deviation*exp(-1) ));
    tau1 = t(idx_t);
    
    for i_t = 1:length(t_a)
        x(i_t,:) = expAt{i_t}*x_init2;
    end
    deviation = max(abs(x(:,2)+f_ss(i)-1));
    
    idx_t = max(find( abs(x(:,2)+f_ss(i)-1) > deviation*exp(-1) ));
    tau2 = t(idx_t);
    
    for i_t = 1:length(t_a)
        x(i_t,:) = expAt{i_t}*x_init3;
    end
    deviation = max(abs(x(:,2)+f_ss(i)-1));
    
    idx_t = max(find( abs(x(:,2)+f_ss(i)-1) > deviation*exp(-1) ));
    tau3 = t(idx_t);
    
    tau(i) = (tau1+tau2+tau3)/3;
end
toc

figure;
plot(V,f_ss);

figure;
plot(V,tau);