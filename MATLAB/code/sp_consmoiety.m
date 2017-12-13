function K = sp_consmoiety(N)

tol = 1e-10;

% Use methods from Schuster and Hofer (1991)
r = size(N,2);

T = [N eye(size(N,1))];

for j = 0:r-1

    S = T(:,r+1:end) == 0;
    T_new = [];
    num_rows = size(T,1);
    for i = 1:num_rows
        for k = i:num_rows
            p = T(i,j+1)*T(k,j+1);
            if p < 0
                Sik = S(i,:).*S(k,:);
                sum_Sik = sum(Sik);
                flag_set = true;
                for l = setdiff(1:num_rows,[i k])
                    Sl = S(l,:);
                    if Sik*Sl' == sum_Sik
                        flag_set = false;
                        break;
                    end
                end
                if flag_set
                    new_vec = abs(T(i,j+1))*T(k,:) + abs(T(k,j+1))*T(i,:);
                    T_new = [T_new; new_vec];
                end
            end
        end
        if abs(T(i,j+1)) < tol
            T(i,j+1) = 0;
            T_new = [T_new; T(i,:)];
        end
    end
    
    T = T_new;

end

K = T(:,r+1:end);

for i = 1:size(K,1)
    K(i,:) = K(i,:)/gcd_array(K(i,:));
end


% Apply modifications in Scuster and Hilgetag (1995)
% num_moieties = size(K,1);
% for i = 1:num_moieties
%     K_i_zero = K(i,:) == 0;
%     for k = i:num_moieties
%         K_k_zero = K(k,:) == 0
%         if any((K_i_zero + K_k_zero) == 1)
%             i
%             k
%         end
%     end
% end

[~,rowind] = rref(K');
K = K(rowind,:);

end
