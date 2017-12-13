function dxdt = df_BG(x,t,t_vec,V_vec,params)
R = 8.314;
T = 310;
F = 96485;

V = interp1(t_vec,V_vec,t);

alpha1 = params(1)*exp(params(2)*F*V/R/T);
beta1 = params(3)*exp(params(4)*F*V/R/T);
alpha2 = params(5)*exp(params(6)*F*V/R/T);
beta2 = params(7)*exp(params(8)*F*V/R/T);
k3f = beta1*alpha2/alpha1/beta2 * 1e5;
k3r = 1e5;

v1 = alpha1*x(2) - beta1*x(1);
v2 = alpha2*x(2) - beta2*x(3);
v3 = k3f*x(1) - k3r*x(3);

dxdt = [v1-v3; -v1-v2; v2+v3];

end