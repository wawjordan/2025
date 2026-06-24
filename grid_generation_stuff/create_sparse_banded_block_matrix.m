function A = create_sparse_banded_block_matrix(m,n,b,data)
D = concatenate_data(m,n,b,data);
[I,J] = all_blk_diag_idx(m,n,b);
A = sparse(I,J,D);
end

function D = concatenate_data(m,n,b,data)
nnz = n^2 *(m+2*(sum(m-(1:(b-1)/2))));
if isscalar(data)
    if isnumeric(data)
        D = data*ones(nnz,1);
        return
    elseif iscell(data)
        if isnumeric(data{1})
            D = [data{:}];
            return
        else
            error('data should either be scalar or cell array of numeric arrays')
        end
    end
end

D = cell(b,1);
for k = 1:b
    d = k-(b+1)/2;
    n_blk_rows     = m-abs(d);
    if any( size(data{k}) ~= [n,n,n_blk_rows])
        error('input size for data{%d}',k)
    end
    D{k} = data{k}(:).';
end
D = [D{:}];
end

function [I,J] = all_blk_diag_idx(m,n,b)
I = cell(b,1);
J = cell(b,1);
for k = 1:b
    d = k-(b+1)/2;
    [I{k},J{k}] = blk_diag_idx(m,n,d);
end
I = [I{:}];
J = [J{:}];
end

function [I,J] = blk_diag_idx(m,n,d)
    n_blk_rows     = m-abs(d);
    I = zeros(n,n,n_blk_rows);
    J = zeros(n,n,n_blk_rows);
    k1 = max(1,1-d);
    k2 = min(m,m-d);
    % for k = max(1,1-d):min(m,m-d)
    for k = k1:k2
        i1 = n*(k-1)+1;
        i2 = i1 + n - 1;
        j1 = n*(k+d-1)+1;
        j2 = j1 + n - 1;
        [I(:,:,k+1-k1),J(:,:,k+1-k1)] = ndgrid(i1:i2,j1:j2);
    end
    I = I(:).';
    J = J(:).';
end