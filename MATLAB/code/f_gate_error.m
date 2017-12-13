function result = f_gate_error(params,V,f_ss,tau_f,compare)
R = 8.314;
T = 310;
F = 96485;

f_ss_fit = f_gate_func(params,V/1000);

error_ss = square_error(log(f_ss_fit) - log(f_ss));

tspan = [0 100];
t = 0:1:100;
t2 = 0:1:200;

f_ss1 = f_ss(V == -80);
tau_f1 = tau_f(V == -80);
fode1 = @(t,x) 1/tau_f1 * (f_ss1 - x);
f_01 = 0;
sol1 = ode15s(fode1,tspan,f_01);
f1 = deval(sol1,t);

f_ss2 = f_ss(V == 40);
tau_f2 = tau_f(V == 40);
fode2 = @(t,x) 1/tau_f2 * (f_ss2 - x);
f_02 = 0;
% sol2a = ode15s(fode2,2*tspan,f_02);
% f2a = deval(sol2a,t2);
sol2b = ode15s(fode2,tspan,1-f_02);
f2b = deval(sol2b,t);

f_ss3 = f_ss(V == 0);
tau_f3 = tau_f(V == 0);
fode3 = @(t,x) 1/tau_f3 * (f_ss3 - x);
f_03 = 1;
sol3 = ode15s(fode3,tspan,f_03);
f3 = deval(sol3,t);

f_ss4 = f_ss(V == -40);
tau_f4 = tau_f(V == -40);
fode4 = @(t,x) 1/tau_f4 * (f_ss4 - x);
f_04a = 1;
f_04b = 0;
sol4a = ode15s(fode4,tspan,f_04a);
f4a = deval(sol4a,t);
sol4b = ode15s(fode4,tspan,f_04b);
f4b = deval(sol4b,t);

alpha1_0 = params(1);
beta1_0 = params(3);

z_f1 = params(2);
z_r1 = params(4);


alpha2_0 = params(5);
beta2_0 = params(7);

z_f2 = params(6);
z_r2 = params(8);

alpha1 = alpha1_0*exp(z_f1*-0.08*F/R/T);
beta1 = beta1_0*exp(z_r1*-0.08*F/R/T);

alpha2 = alpha2_0*exp(z_f2*-0.08*F/R/T);
beta2 = beta2_0*exp(z_r2*-0.08*F/R/T);
sim_ode1 = @(t,x) dx_f2(x,alpha1,beta1,alpha2,beta2);
x_01 = [0; 1; 0];

alpha1 = alpha1_0*exp(z_f1*0.04*F/R/T);
beta1 = beta1_0*exp(z_r1*0.04*F/R/T);

alpha2 = alpha2_0*exp(z_f2*0.04*F/R/T);
beta2 = beta2_0*exp(z_r2*0.04*F/R/T);
sim_ode2 = @(t,x) dx_f2(x,alpha1,beta1,alpha2,beta2);
x_02a = [0; 1; 0];
x_02b = [0; 0; 1];

alpha1 = alpha1_0;
beta1 = beta1_0;

alpha2 = alpha2_0;
beta2 = beta2_0;
sim_ode3 = @(t,x) dx_f2(x,alpha1,beta1,alpha2,beta2);
x_03a = [1; 0; 0];
x_03b = [0; 0; 1];

alpha1 = alpha1_0*exp(z_f1*-0.04*F/R/T);
beta1 = beta1_0*exp(z_r1*-0.04*F/R/T);

alpha2 = alpha2_0*exp(z_f2*-0.04*F/R/T);
beta2 = beta2_0*exp(z_r2*-0.04*F/R/T);
sim_ode4 = @(t,x) dx_f2(x,alpha1,beta1,alpha2,beta2);
x_04a = [1; 0; 0];
x_04b = [0; 1; 0];

sol1 = ode15s(sim_ode1,tspan,x_01);
x1 = deval(sol1,t,2);

% sol2a = ode15s(sim_ode2,2*tspan,x_02a);
% x2a = deval(sol2a,t2,2);

sol2b = ode15s(sim_ode2,tspan,x_02b);
x2b = deval(sol2b,t,2);

sol3a = ode15s(sim_ode3,tspan,x_03a);
x3a = deval(sol3a,t,2);

sol3b = ode15s(sim_ode3,tspan,x_03b);
x3b = deval(sol3b,t,2);

sol4a = ode15s(sim_ode4,tspan,x_04a);
x4a = deval(sol4a,t,2);

sol4b = ode15s(sim_ode4,tspan,x_04b);
x4b = deval(sol4b,t,2);

error_1 = square_error((1-x1)-f1);
% error_2a = square_error((1-x2a)-f2a);
error_2b = square_error((1-x2b)-f2b);
error_3a = square_error(log(1-x3a)-log(f3));
error_3b = square_error(log(1-x3b)-log(f3));
error_4a = square_error((1-x4a)-f4a);
error_4b = square_error((1-x4b)-f4b);

if compare
    figure;
    plot(t,f1,t,1-x1,'LineWidth',2);
    xlabel('Time (ms)');
    ylabel('f');
    set(gca,'FontSize',16);

    figure;
    plot(t,f2b,t,1-x2b,'LineWidth',2);
    xlabel('Time (ms)');
    ylabel('f');
    set(gca,'FontSize',16);
    
    figure;
    plot(t,f3,t,1-x3a,'LineWidth',2);
    xlabel('Time (ms)');
    ylabel('f');
    set(gca,'FontSize',16);
    
    figure;
    plot(t,f3,t,1-x3b,'LineWidth',2);
    xlabel('Time (ms)');
    ylabel('f');
    set(gca,'FontSize',16);
    
    figure;
    plot(t,f4a,t,1-x4a,'LineWidth',2);
    xlabel('Time (ms)');
    ylabel('f');
    set(gca,'FontSize',16);
    
    figure;
    plot(t,f4b,t,1-x4b,'LineWidth',2);
    xlabel('Time (ms)');
    ylabel('f');
    set(gca,'FontSize',16);
end

result = 5*error_ss + error_1 + 3*error_2b + error_3a + 3*error_3b + error_4a + error_4b;


end