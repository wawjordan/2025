function [pexp,pidx] = pexponents_alt_old(dimen,total_degree)
n_unknowns = nchoosek(dimen+total_degree,total_degree);
pexp = zeros(dimen,n_unknowns);
full_exp = total_degree*ones(1,dimen);
pidx = zeros(n_unknowns,1);
full_degree = (total_degree+1)^dimen;
cnt = 0;
nSub = (total_degree+1)*ones(dimen,1);
p_idx = zeros(n_unknowns,1);
n10s = 10^floor(log10(full_degree));
for n = 0:full_degree
    tmp_exp = global_to_local(n+1,nSub) - 1;
    current_degree = sum(tmp_exp);
    if current_degree <= total_degree
        cnt = cnt + 1;
        pidx(cnt) = local_to_global(tmp_exp+1,full_exp+1);
        pexp(:,cnt) = tmp_exp;
        % (order,reverse lexicographic in each dim) sort
        p_idx(cnt) = current_degree*n10s*10^(dimen-1);
        for j = 1:dimen
            p_idx(cnt) = p_idx(cnt) + (total_degree-pexp(j,cnt))*10^(dimen-j);
        end
    end
end

% [~,idx] = sort(p_idx);
% pexp = pexp(:,idx);
end

function iSub = global_to_local(iG,nSub)
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

function iG = local_to_global(iSub,nSub)
iSub = iSub(:);
nSub = nSub(:);
nDims = numel(iSub);
p = 1;
iG = 1;
for i = 1:nDims
    iG = iG + ( iSub(i) - 1 )*p;
    p = p*nSub(i);
end
end
