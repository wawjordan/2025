function [M] = computational_transform_matrix(X,degree,scale)
flag = false;
if ~isa(X,"cell")
    flag = true;
end
n_dim = numel(X);
if (n_dim>3)
    error('only dimensions 1-3 are supported')
end
for i = 1:n_dim
    if ~isa(X{i},'numeric')
        flag = true;
        break;
    elseif ~isvector(X{i})
        flag = true;
        break;
    end
end
if (flag)
    error('X must be a cell array of 1D monotonic vectors')
end
switch n_dim
    case(1)
        xvec = ( X{1}(:) ).';
    case(2)
        [x1,x2] = ndgrid( X{1}, X{2} );
        xvec = [x1(:),x2(:)].';
    case(3)
        [x1,x2,x3] = ndgrid( X{1}, X{2}, X{3} );
        xvec = [x1(:),x2(:),x3(:)].';
end
xvec = xvec*scale;
% [exponents,~,~] = get_exponents( n_dim, degree );
exponents = get_exponents_tp( n_dim, degree );
n_terms = size(exponents,2);
n_nodes = prod( cellfun(@numel,X) );

M = ones(n_nodes,n_terms);
for j = 1:n_terms
    for i = 1:n_nodes
        M(i,j) = evaluate_monomial( n_dim, exponents,j,xvec(:,i));
    end
end
end

function [val,coef] = evaluate_monomial( n_dim, exponents, term, x )
val = 1.0;
coef = 1;
for d = 1:n_dim
    for i = exponents(d,term):-1:1
        val = val*x(d);
        coef = coef * i;
    end
end
end

function [dval,dcoef,coef] = evaluate_monomial_derivative( n_dim,     ...
                                                           exponents, ...
                                                           term,      ...
                                                           x,         ...
                                                           order )
dval  = 0.0;
dcoef = 1;
coef  = 1;
if ( any( exponents(:,term) - order(:) < 0 ) ); return; end
dval = 1.0;
for d = 1:n_dim
    for i = exponents(d,term):-1:exponents(d,term)-order(d)+1
        dcoef = dcoef * i;
    end
    for i = exponents(d,term)-order(d):-1:1
        coef = coef * i;
        dval = dval*x(d);
    end
end
end

function exponents = get_exponents_tp( n_dim, degree )
n_terms   = (degree+1)^n_dim;
exponents = zeros(n_dim,n_terms);
nsub = (degree+1)*ones(1,n_dim);
for j = 1:n_terms
    exponents(:,j) = global2local(j,nsub)-1;
end
if (n_dim== 2 && degree==2)
    exponents = exponents(:,[1,2,4,5,3,7,9,6,8]);
end
% % 1: [0,0]  1
% % 2: [1,0]  2
% % 3: [0,1]  4
% % 4: [1,1]  5
% % 5: [2,0]  3
% % 6: [0,2]  7
% % 7: [2,2]  9
% % 8: [2,1]  6
% % 9: [1,2]  8
end

% function exponents = get_exponents_tp( n_dim, degree )
% max_degree   = degree*n_dim;
% [exponents,~,~] = get_exponents( n_dim, max_degree );
% exponents(:,any(exponents>degree,1)) = [];
% % 1: [0,0]  1
% % 2: [1,0]  2
% % 3: [0,1]  3
% % 4: [1,1]  5
% % 5: [2,0]  4
% % 6: [0,2]  6
% % 7: [2,2]  9
% % 8: [2,1]  7
% % 9: [1,2]  8
% end

function [exponents,idx,diff_idx] = get_exponents( n_dim, degree )
n_terms   = nchoosek( n_dim + degree, degree );
exponents = zeros(n_dim,n_terms);
idx       = zeros(degree+1,1);
diff_idx  = zeros(n_dim,n_terms);
cnt = 0;
for curr_total_degree = 0:degree
    nsub(1:n_dim) = curr_total_degree + 1;
    N_full_terms = (curr_total_degree+1)^n_dim;
    for j = 0:N_full_terms
        tmp_exp = global2local(j+1,nsub)-1;
        if ( sum(tmp_exp) == curr_total_degree )
            cnt = cnt + 1;
            exponents(:,cnt) = tmp_exp;
        end
    end
    idx(curr_total_degree+1) = cnt;
end
diff_idx(1:n_dim,1:n_terms) = -1;
if ( degree == 0 ); return; end
for j = 1:idx(degree)
    tmp_exp = exponents(:,j);
    curr_total_degree = sum(tmp_exp);
    cnt = 0;
    for i = idx(curr_total_degree+1)+1:idx(curr_total_degree+2)
        if ( sum( abs(exponents(:,i) - tmp_exp) ) == 1 )
            cnt = cnt +1;
            diff_idx(cnt,j) = i;
        end
    end
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

% function iG = local2global(iSub,nSub)
% iSub = iSub(:);
% nSub = nSub(:);
% nDims = numel(iSub);
% p = 1;
% iG = 1;
% for i = 1:nDims
%     iG = iG + ( iSub(i) - 1 )*p;
%     p = p*nSub(i);
% end
% end