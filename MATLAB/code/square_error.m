function result = square_error(vector)
    if ~iscolumn(vector)
        vector = transpose(vector);
    end
    result = transpose(vector)*vector;
end