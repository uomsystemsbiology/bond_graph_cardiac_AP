function peak_Ca = calc_peak_Ca(t,Cai)

num_beats = floor(t(end));
APD = zeros(num_beats,1);

for i_beat = 1:num_beats
    Idx_beat = ((1+1000*(i_beat-1)):(1+1000*i_beat));
    Ca_beat = Cai(Idx_beat);
    peak_Ca(i_beat) = max(Ca_beat);
end



end