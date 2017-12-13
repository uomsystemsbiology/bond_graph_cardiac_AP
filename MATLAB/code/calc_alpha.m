function alpha = calc_alpha(params,V)
% params: [kf, zf, kr, zr];
R = 8.314;
T = 310;
F = 96485;

alpha = params(1)*exp(params(2)*F*V/R/T);

end