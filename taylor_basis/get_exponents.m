function [exp,idx,exp_show] = get_exponents( n_dim, degree )
n_terms = get_n_terms(n_dim,degree);
exp     = zeros(n_dim,n_terms);
idx     = zeros(degree+1,1);
cnt = 0;
for curr_total_degree = 0:degree
    N_full_terms = (curr_total_degree+1)^n_dim;
    for j = 0:N_full_terms
        nsub = (curr_total_degree + 1)*ones(1,n_dim);
        tmp_exp = global2local(j+1,nsub)-1;
        if ( sum(tmp_exp) == curr_total_degree )
            cnt = cnt + 1;
            exp(:,cnt) = tmp_exp;
        end
    end
    idx(curr_total_degree+1) = cnt;
end
exp_show = [(1:n_terms).',exp.'];
end

function n_terms = get_n_terms( dimen, total_degree )
if (total_degree < 0)
    n_terms = 0;
else
    n_terms = nchoosek( dimen + total_degree, total_degree );
end
end

function iSub = global2local(iG,nSub)
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