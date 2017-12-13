function T = calcT(I_vec,num_rows)
num_cols = length(I_vec);
T = zeros(num_rows,num_cols);
for i = 1:num_cols
    T(I_vec(i),i) = 1;
end
end