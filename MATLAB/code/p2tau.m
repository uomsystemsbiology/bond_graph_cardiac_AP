function tau = p2tau(params,V)
[alpha,beta] = calc_gate_transitions(params,V/1000);
tau = calc_tau(alpha,beta);
end