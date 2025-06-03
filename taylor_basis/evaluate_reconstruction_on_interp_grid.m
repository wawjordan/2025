function [X,F] = evaluate_reconstruction_on_interp_grid(coefs,T,Q,npts)

% get dimension
n_dim = T.n_dim;

X = get_local_interp_grid(Q,n_dim,npts);
siz = num2cell(size(X{1}));
siz2 = cellfun(@(x)ones(x,1),siz,'UniformOutput',false);
siz2 = [siz2,3];
Xtmp = cat(n_dim+1,X{:});
Xtmp = mat2cell(Xtmp,siz2{:});
    
coefs_in = zeros(T.n_terms,1);
coefs_in(1:numel(coefs)) = coefs;

F = cellfun(@(X)T.reconstruct_from_coefs(coefs_in,[X(:)]),Xtmp);

end