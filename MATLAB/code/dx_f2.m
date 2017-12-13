function dx = dx_f2(x,alpha1,beta1,alpha2,beta2)

v1 = alpha1*x(2) - beta1*x(1);
v2 = alpha2*x(2) - beta2*x(3);

dx = [v1; -(v1+v2); v2];
end