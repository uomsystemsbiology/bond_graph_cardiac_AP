function beta = calc_beta(params,V)
% params: [kf, zf, kr, zr];
R = 8.314;
T = 310;
F = 96485;

beta = params(3)*exp(params(4)*F*V/R/T);

end