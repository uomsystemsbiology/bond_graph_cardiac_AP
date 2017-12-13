function [Nf,Nr] = load_stoichiometry(data_dir,sys_name)
forward_mat_path = [data_dir sys_name '_forward_matrix.txt'];
reverse_mat_path = [data_dir sys_name '_reverse_matrix.txt'];
Nf = load_matrix(forward_mat_path);
Nr = load_matrix(reverse_mat_path);
end