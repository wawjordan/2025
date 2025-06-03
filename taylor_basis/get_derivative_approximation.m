function dvals = get_derivative_approximation(T_array,vals)
n_nbors = numel(T_array);
n_dim = T_array(1).n_dim;
M = zeros(n_nbors-1,n_dim);
R = zeros(n_nbors-1,1);
for j = 2:n_nbors
    M(j-1,:) = T_array(1).build_linear_reconstruction_matrix( T_array(j) );
    % [T_array(1).moments,T_array(j).moments]
    R(j-1)   = T_array(1).build_linear_reconstruction_rhs(    T_array(j), vals(1), vals(j) );
end
% M
% R
dvals = M\R;

end