function vss = NCX_vss_fitting(params_vec,struct_input,e_0)
R = 8.314;
T = 310;
F = 96485;

W_i = 38;
W_e = 5.182;

kappa_3 = params_vec(1);
kappa_6 = params_vec(2);
K_1 = params_vec(3);
K_2 = params_vec(4);
K_3 = params_vec(5);
K_4 = params_vec(6);
K_5 = params_vec(7);
K_6 = params_vec(8);

K_Nai = params_vec(9);
K_Cai = params_vec(10);
K_Nae = W_i*K_Nai/W_e;
K_Cae = W_i*K_Cai/W_e;

Delta = params_vec(11);

V = struct_input.V;

gam_Cae = K_Cae*W_e*struct_input.Cae;
gam_Cai = K_Cai*W_i*struct_input.Cai;
gam_Nae = K_Nae*W_e*struct_input.Nae;
gam_Nai = K_Nai*W_i*struct_input.Nai;
gam_mem = exp(F*V/R/T);

vss = NCX_vss(Delta,K_1,K_2,K_3,K_4,K_5,K_6,e_0,gam_Cae,gam_Cai,gam_Nae,gam_Nai,gam_mem,kappa_3,kappa_6);

end