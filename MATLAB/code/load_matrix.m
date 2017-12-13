function matrix = load_matrix(stoich_path)

stoich_file_id = fopen(stoich_path,'r');

stoich_file_data = textscan(stoich_file_id,'%s','delimiter','\n');
fclose(stoich_file_id);

num_rows = length(stoich_file_data{1});
num_cols = sum(stoich_file_data{1}{1} == ',')+1;

matrix = zeros(num_rows,num_cols);

for i_row = 1:num_rows
    line_str = stoich_file_data{1}{i_row};
    line_split = regexp(line_str,',','split');
    
    for i_col = 1:num_cols
        matrix(i_row,i_col) = str2double(line_split{i_col});
    end
end


end