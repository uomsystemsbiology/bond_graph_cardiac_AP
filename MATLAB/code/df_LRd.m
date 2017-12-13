function dfdt = df_LRd(f,t,t_vec,V_vec)

V = interp1(t_vec,V_vec,t);

f_ss = 1./(1+exp((V+32)/8)) + 0.6./(1+exp((50-V)/20));
tau_f = 1/(0.0197*exp(-(0.0337*(V+10))^2)+0.02); % Unit ms

dfdt = (f_ss-f)/tau_f;

end