function struct_output = extract_results(path_pref,num_segments,G_charge)

flag_charge = exist('G_charge','var');

C_m = 153400;
W_i = 38;
W_e = 5.182;

t = [];
V = [];
Ki = [];
Cai = [];
Nai = [];
Ke = [];
Cae = [];
Nae = [];
charge = [];

for i_segment = 1:num_segments
    load([path_pref num2str(i_segment) '.mat']);
    t = [t; VOI];
    V = [V; STATES(:,7)/C_m];
    Ki = [Ki; STATES(:,4)/W_i];
    Cai = [Cai; STATES(:,2)/W_i];
    Nai = [Nai; STATES(:,6)/W_i];
    Ke = [Ke; STATES(:,3)/W_e];
    Cae = [Cae; STATES(:,1)/W_e];
    Nae = [Nae; STATES(:,5)/W_e];
    
    if flag_charge
        charge = [charge; double(STATES)*transpose(G_charge)];
    end
end

[~,idx_unique] = unique(t);
struct_output.t = t(idx_unique);
struct_output.V = V(idx_unique);
struct_output.Ki = Ki(idx_unique);
struct_output.Cai = Cai(idx_unique);
struct_output.Nai = Nai(idx_unique);
struct_output.Ke = Ke(idx_unique);
struct_output.Cae = Cae(idx_unique);
struct_output.Nae = Nae(idx_unique);
if flag_charge
    struct_output.charge = charge(idx_unique);
end



end