function [X,F,delta] = evaluate_shape_function(T,N,Q,npts)

% get dimension
n_dim = T.n_dim;

X = get_local_interp_grid(Q,n_dim,npts);

[F,delta] = arrayfun(@(x1,y1)T.eval(N,[x1;y1]),X{1},X{2});

end