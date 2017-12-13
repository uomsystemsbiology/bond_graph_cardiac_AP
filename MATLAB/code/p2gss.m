function g_ss = p2gss(params,V)
[alpha,beta] = calc_gate_transitions(params,V/1000);
g_ss = calc_gss(alpha,beta);
end