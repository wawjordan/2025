%% more grevlex testing (08/22/2025)
clc; clear; close all;


n_dim = 3; degree = 3;

% Nm1 = 0;
% for i = 0:degree
%     N = nchoosek(n_dim+i,i);
%     for j = Nm1:N-1
%         fprintf('%d  ',j+1)
%     end
%     fprintf('\n')
%     Nm1 = N;
% end



% prev_total_degree = 0;
% curr_degree       = 0;
% 
% cnt = 0;
% 
% N_terms_m1 = 0;
% for curr_total_degree = 0:degree
%     N_terms = nchoosek(n_dim+curr_total_degree,curr_total_degree);
%     N_full_terms = (curr_total_degree+1)^n_dim;
%     nsub = repmat(curr_total_degree+1,[1,n_dim]);
%     for j = 0:N_full_terms
%         tmp_exp = zero_mean_basis6.global_to_local(j+1,nsub)-1;
%         curr_degree = sum(tmp_exp);
%         if curr_degree == curr_total_degree
%             cnt = cnt + 1;
%             fprintf('%d  (%d %d %d)\n',cnt,tmp_exp)
%         end
%     end
%     % if cnt == N_terms
%         fprintf('\n')
%     %     N_terms_m1 = N_terms;
%     % end
% end


E = enumerate_exponents(n_dim,degree);
% E.'


E = enumerate_exponents_sorted(n_dim,degree);


function E = enumerate_exponents(n_dim,degree)

E = nan(n_dim, nchoosek(n_dim+degree,degree) );

cnt = 0;
for curr_total_degree = 0:degree
    
    % nsub = repmat(curr_total_degree+1,[1,n_dim]);

    N_full_terms = (curr_total_degree+1)^n_dim;
    for j = 0:N_full_terms
        % tmp_exp1 = zero_mean_basis6.global_to_local(j+1,nsub)-1;
        tmp_exp = nan(n_dim,1);
        p = N_full_terms;
        j_tmp = j+1;
        for d = n_dim:-1:1
            p = fix( p/(curr_total_degree+1) );
            d_tmp = mod(j_tmp-1,p)+1;
            tmp_exp(d) = fix( (j_tmp-d_tmp)/p );
            j_tmp = d_tmp;
        end
        curr_degree = sum(tmp_exp);
        if curr_degree == curr_total_degree
            cnt = cnt + 1;
            % fprintf('%d  (%d %d %d)\n',cnt,tmp_exp)
            E(:,cnt) = tmp_exp;
        end
    end
end

end


function E = enumerate_exponents_sorted(n_dim,degree)

E = nan(n_dim, nchoosek(n_dim+degree,degree) );
% idea was to try and use a swap sort on the array while building it, but that might not be possible
cnt = 0;
for curr_total_degree = 0:degree
    N_full_terms = (curr_total_degree+1)^n_dim;
    for j = 0:N_full_terms
        % tmp_exp1 = zero_mean_basis6.global_to_local(j+1,nsub)-1;
        tmp_exp = nan(n_dim,1);
        p = N_full_terms;
        j_tmp = j+1;
        for d = n_dim:-1:1
            p = fix( p/(curr_total_degree+1) );
            d_tmp = mod(j_tmp-1,p)+1;
            tmp_exp(d) = fix( (j_tmp-d_tmp)/p );
            j_tmp = d_tmp;
        end
        curr_degree = sum(tmp_exp);
        if curr_degree == curr_total_degree
            cnt = cnt + 1;
            E(:,cnt) = tmp_exp;
            if (cnt > 1)
                for d = 1:n_dim
                    if tmp_exp(d) < E(d,cnt-1)
                        E(:,cnt) = E(:,cnt-1);
                        E(:,cnt-1) = tmp_exp;
                    end
                end
            end
        end
    end
end

end