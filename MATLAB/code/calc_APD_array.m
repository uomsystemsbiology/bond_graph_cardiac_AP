function APD = calc_APD_array(t,V)

num_beats = floor(t(end));
APD = zeros(num_beats,1);

for i_beat = 1:num_beats
    Idx_beat = ((1+1000*(i_beat-1)):(1+1000*i_beat));
    t_beat = 1000*(t(Idx_beat)-i_beat+1)-300;
    V_beat = V(Idx_beat);
    V_peak = max(V_beat);
    V_end = V_beat(end);
    V_90 = 0.9*V_end + 0.1*V_peak;
    Idx_repolarised = find((V_beat <= V_90) & (t_beat > 20),1);
    Idx_cross = [Idx_repolarised-1 Idx_repolarised];
    APD(i_beat) = interp1(V_beat(Idx_cross),t_beat(Idx_cross),V_90);
end



end