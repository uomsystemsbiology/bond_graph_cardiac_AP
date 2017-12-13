function I_GHK = calc_IGHK(G_GHK,V,cKi_st,cKo)
R = 8.314;
T = 310;
F = 96485;

V_norm = F*V/R/T;
GHK_factor = V_norm./(1-exp(-V_norm));
idx_nan = isnan(GHK_factor);
GHK_factor(idx_nan) = 1;
I_GHK = G_GHK*GHK_factor.*(cKi_st - cKo*exp(-V_norm));

end