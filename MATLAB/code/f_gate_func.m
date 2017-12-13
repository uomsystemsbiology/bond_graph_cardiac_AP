function [f_ss,tau] = f_gate_func(params,V)
R = 8.314;
T = 310;
F = 96485;

alpha1_0 = params(1);
beta1_0 = params(3);

z_f1 = params(2);
z_r1 = params(4);


alpha2_0 = params(5);
beta2_0 = params(7);

z_f2 = params(6);
z_r2 = params(8);


x_init1 = [1; 0; 0];
x_init2 = [0; 1; 0];
x_init3 = [0; 0; 1];

options = odeset('NonNegative',1:3,'RelTol',1e-3,'AbsTol',1e-8,'Refine',100);
tspan = [0 100];

f_ss = zeros(length(V),1);
tau = zeros(length(V),1);

for i = 1:length(V)
    alpha1 = alpha1_0*exp(z_f1*V(i)*F/R/T);
    beta1 = beta1_0*exp(z_r1*V(i)*F/R/T);
    K1 = alpha1/beta1;

    alpha2 = alpha2_0*exp(z_f2*V(i)*F/R/T);
    beta2 = beta2_0*exp(z_r2*V(i)*F/R/T);
    K2 = alpha2/beta2;
    
    f_ss(i) = (K1+K2)/(1+K1+K2);
    
    if nargout > 1
        sim_ode = @(t,x) dx_f2(x,alpha1,beta1,alpha2,beta2);
        [t,x] = ode15s(sim_ode,tspan,x_init1,options);
        deviation = max(abs(x(:,2)+f_ss(i)-1));

        idx_t = find( abs(x(:,2)+f_ss(i)-1) > deviation*exp(-1) ,1,'last');
        tau1 = t(idx_t);
        if isempty(tau1)
            tau1 = tspan(2);
        end

        [t,x] = ode15s(sim_ode,tspan,x_init2,options);
        deviation = max(abs(x(:,2)+f_ss(i)-1));

        idx_t = find( abs(x(:,2)+f_ss(i)-1) > deviation*exp(-1) ,1,'last');
        tau2 = t(idx_t);
        if isempty(tau2)
            tau2 = tspan(2);
        end

        [t,x] = ode15s(sim_ode,tspan,x_init3,options);
        deviation = max(abs(x(:,2)+f_ss(i)-1));

        idx_t = find( abs(x(:,2)+f_ss(i)-1) > deviation*exp(-1) ,1,'last');
        tau3 = t(idx_t);
        if isempty(tau3)
            tau3 = tspan(2);
        end

        tau(i) = (tau1+tau2+tau3)/3;
    end
end