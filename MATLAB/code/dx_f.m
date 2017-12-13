function dx = dx_f(x,V)

alpha1_0 = 5;
beta1_0 = exp(3);

z_f1 = 10;
z_r1 = 100;

alpha2_0 = 5;
beta2_0 = exp(3);

z_f2 = -10;
z_r2 = -100;

alpha1 = alpha1_0*exp(z_f1*V);
beta1 = beta1_0*exp(z_r1*V);

alpha2 = alpha2_0*exp(z_f2*V);
beta2 = beta2_0*exp(z_r2*V);

v1 = alpha1*x(2) - beta1*x(1);
v2 = alpha2*x(2) - beta2*x(3);

dx = [v1; -(v1+v2); v2];
end