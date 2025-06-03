function [X,F] = evaluate_function_on_interp_grid(fun,Q,npts)
n_dim = nargin(fun);
X = get_local_interp_grid(Q,n_dim,npts);

F = fun(X{1:n_dim});
end