function total_error = error_BW(params_vec,struct_data)
% Experimental conditions:
struct_input = struct('Nai',15,...
    'Cae',2,...
    'Cai',450e-6,...
    'V',struct_data.V(end),...
    'Nae',135);

% Normalising factors
normalising_factor_data = abs(struct_data.array_I(end));
normalising_factor_model = -abs(NCX_vss_fitting(params_vec,struct_input,1));

total_error = 0;
for i_V = 1:length(struct_data.V)
    struct_input.V = struct_data.V(i_V);
    r_I_model = NCX_vss_fitting(params_vec,struct_input,1)/normalising_factor_model;
    r_I_data = struct_data.array_I(i_V)/normalising_factor_data;
    if struct_data.V(i_V) >= -0.12
        total_error = total_error + (r_I_model-r_I_data)^2;
    end
end

end