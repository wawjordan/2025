function avg = get_cell_avg_of_function(Q,fun)
avg = arrayfun(@(Q)get_single_cell_avg(Q,fun,nargin(fun)),Q);
end

function avg = get_single_cell_avg(Q,fun,n_dim)
tmp_val = zeros(1,Q.n_quad);
for n = 1:Q.n_quad
    tmp_x      = num2cell(Q.quad_pts(1:n_dim,n));
    tmp_val(n) = fun(tmp_x{:});
end
volume = Q.integrate(ones(1,Q.n_quad));
avg    = Q.integrate(tmp_val)/volume;
end