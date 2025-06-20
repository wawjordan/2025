function moment_hat = compute_shifted_moment_mod_old( pexp, xi, xj, exponents, moments, p )
n_dim = numel(pexp);
deltas = zeros(n_dim,p.Nmax);
deltas(:,1) = 1.0;
for m = 1:max(pexp)
    deltas(:,m+1) = (xj(1:n_dim) - xi(1:n_dim)).^m;
end

moment_hat = 0;
for m = 1:prod(pexp+1)
    cexp = global_to_local_local(m,pexp+1)-1;
    
    comb_factor = 1.0;
    for j = 1:n_dim
        comb_factor = comb_factor*p.binomial_coefficient(pexp(j),cexp(j));
    end
    grid_factor = 1.0;
    for j = 1:n_dim
        grid_factor = grid_factor*deltas(1,cexp(j)+1);
    end

    % n = local_to_global(pexp-cexp+1);
    % for shift_row = 1:numel(pidx)
    %     if ( n == pidx(shift_row) )
    %         break
    %     end
    % end
    shift_exp = pexp-cexp;
    for shift_row = 1:numel(moments)
        if all( shift_exp == exponents(:,shift_row) )
            break
        end
    end
    shift_moment = moments(shift_row);

    moment_hat = moment_hat + comb_factor * grid_factor * shift_moment;
end

end

function iSub = global_to_local_local(iG,nSub)
nDims = numel(nSub);
iSub = zeros(1,nDims);
if (nDims==1)
    iSub(1) = iG;
    return
end
p = prod(nSub);
iGtmp = iG;
for i = nDims:-1:1
    p = fix( p/nSub(i) );
    iTmp = mod(iGtmp-1,p)+1;
    iSub(i) = fix( (iGtmp-iTmp)/p ) + 1;
    iGtmp = iTmp;
end

end