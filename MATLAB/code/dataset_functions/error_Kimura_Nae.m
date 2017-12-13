function total_error = error_Kimura_Nae(params_vec,struct_data)
% Experimental conditions:
struct_input = struct('Nai',0,...
    'Cae',1,...
    'Cai',430e-6,...
    'V',struct_data.V(1),...
    'Nae',struct_data.Nae(4));

% Normalising factors
normalising_factor_data = abs(struct_data.array_I(4,1));
normalising_factor_model = -abs(NCX_vss_fitting(params_vec,struct_input,1));

total_error = 0;
for i_Nae = 1:length(struct_data.Nae)
    struct_input.Nae = struct_data.Nae(i_Nae);
    for i_V = 1:length(struct_data.V)
        struct_input.V = struct_data.V(i_V);
        r_I_model = NCX_vss_fitting(params_vec,struct_input,1)/normalising_factor_model;
        r_I_data = struct_data.array_I(i_Nae,i_V)/normalising_factor_data;
        total_error = total_error + (r_I_model-r_I_data)^2;
    end
end

end