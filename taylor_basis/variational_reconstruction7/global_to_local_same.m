function iSub = global_to_local_same(iG,n_dim,n_max)

% nDims = numel(nSub);
iSub = zeros(1,n_dim);
if (n_dim==1)
    iSub(1) = iG;
    return
end
% p = prod(nSub);
p = n_max^n_dim;
iGtmp = iG;
for i = n_dim:-1:1
    p = fix( p/n_max );
    iTmp = mod(iGtmp-1,p)+1;
    iSub(i) = fix( (iGtmp-iTmp)/p ) + 1;
    iGtmp = iTmp;
end
end