%% Polynomial composed with affine map (09/11/2026)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
n_dim  = 2;
degree = 5;
term = 18;
[exponents,idx,diff_idx] = get_exponents( n_dim, degree );

syms x b [n_dim 1] real
syms A [n_dim n_dim] real

px = prod(x.^exponents(:,term),1);
% px = sum(prod(x.^exponents,1));
y = A*x + b;
py = subs(px,sym2cell(x),sym2cell(y));

g = gradient(py,x);

dp1_1 = diff_poly(py,x,[1,0]);

str = clean_fortran_sym(simplify(dp1_1)-1e-2);
% dp1_2 = diff_poly(py,x,[0,1]);

dp2_1 = diff_affine_poly(px,x,A,b,[1,0]);
dp3_1 = diff_affine_poly2(px,x,A,b,[1,0]);

% dp2_2 = diff_affine_poly(px,x,A,b,[0,1]);

for j = 0:degree
    for i = 0:degree
        dp1 = diff_poly(py,x,[i,j]);
        dp2 = diff_affine_poly(px,x,A,b,[i,j]);
        % dp3 = diff_affine_poly2(px,x,A,b,[i,j])
        check = logical( simplify(dp1==dp2) );
        if ~check
            error('derivative does not match')
        end
    end
end



% T1 = tensor_diff(px,x,1); check_tensor_diff(T1,px,x);
% T2 = tensor_diff(px,x,2); check_tensor_diff(T2,px,x);
% all(logical(simplify(hessian(px,x)==T2)))
% T3 = tensor_diff(px,x,3); check_tensor_diff(T3,px,x);
% T4 = tensor_diff(px,x,4); check_tensor_diff(T4,px,x);
% T5 = tensor_diff(px,x,5); check_tensor_diff(T5,px,x);

function dp = diff_poly(p,x,order)
n_dim = numel(x);
dp = p;
for d = 1:n_dim
    dp = diff(dp,x(d),order(d));
end
end

function dp = diff_affine_poly(p,x,A,b,order)
n_dim = numel(x);
dp = p;
for d = 1:n_dim
    col = d;
    for o = 1:order(d)
        u = A(:,col);
        dp = directional_derivative(dp,x,u);
    end
end
y = A*x + b;
dp = subs(dp,sym2cell(x),sym2cell(y));
function dp = directional_derivative(p,x,u)
dp = dot(gradient(p,x),u);
end
end

function dp = diff_affine_poly2(p,x,A,b,order)
n_dim = numel(x);


total_order = sum(order);
idx1 = order2idx(order);
if (total_order==0)
    dp = p;
    n_terms = 0;
else
    dp = sym(0);
    n_terms = n_dim^total_order;
end

nsub = n_dim*ones(1,total_order);


for i = 1:n_terms
    idx = global2local(i,nsub);
    xs  = sym2cell(x(idx));
    tmp = diff(p,xs{:});
    if logical(tmp==0)
        continue
    end
    % product of A terms
    Atmp = sym(1);
    for k = 1:total_order
        ii = idx(k);
        jj = idx1(k);
        Atmp = Atmp*A(ii,jj);
    end
    dp = dp + Atmp*tmp;
end
% y = A*x + b;
% dp = subs(dp,sym2cell(x),sym2cell(y));
end

function idx = order2idx(order)
% total order == 5
% [5,0] --> [1,1,1,1,1] --> 1  = 1
% [4,1] --> [1,1,1,1,2] --> 17 = 1 + 2^4
% [3,2] --> [1,1,1,2,2] --> 25 = 1 + 2^4 + 2^3
% [2,3] --> [1,1,2,2,2] --> 29 = 1 + 2^4 + 2^3 + 2^2
% [1,4] --> [1,2,2,2,2] --> 31 = 1 + 2^4 + 2^3 + 2^2 + 2^1
% [0,5] --> [2,2,2,2,2] --> 32 = 1 + 2^4 + 2^3 + 2^2 + 2^1 + 2^0
total_order = sum(order);
idx = zeros(1,total_order);
cnt = 0;
for d = 1:numel(order)
    for o = 1:order(d)
        cnt = cnt + 1;
        idx(cnt) = d;
    end
end
end


function dp = diff_affine_poly3(p,x,A,b,order)
n_dim = numel(x);


total_order = sum(order);
if (total_order==0)
    dp = p;
else
    dp = sym(0);
end
n_terms = n_dim^total_order;
nsub = n_dim*ones(1,total_order);
tmp = sym(zeros(n_terms,1));
for i = 1:n_terms
    idx = global2local(i,nsub);
    xs  = sym2cell(x(idx));
    tmp(i) = diff(p,xs{:});
end


for i = 1:n_terms
    idx = global2local(i,nsub);
    for d = 1:n_dim
        col = d;
        for o = 1:order(d)
            % product of A terms
            % Atmp = tmp(i);
            Atmp = sym(1);
            for k = 1:total_order
                ii = idx(k);
                Atmp = Atmp*A(ii,col);
            end
            dp = dp + Atmp*tmp(i);
        end
    end
end

% for d = 1:n_dim
%     col = d;
%     for o = 1:order(d)
%         for i = 1:n_terms
%             idx = global2local(i,nsub);
%             % product of A terms
%             % Atmp = tmp(i);
%             Atmp = sym(1);
%             for k = 1:total_order
%                 ii = idx(k);
%                 Atmp = Atmp*A(ii,col);
%             end
%             dp = dp + Atmp*tmp(i);
%         end
%     end
% end

% for j = 1:total_order
%     for i = 1:n_terms
%         idx = global2local(i,nsub);
%         % product of A terms
%         % Atmp = tmp(i);
%         Atmp = sym(1);
%         for k = 1:total_order
%             ii = idx(k);
%             Atmp = Atmp*A(ii,j);
%         end
%         dp = dp + Atmp*tmp(i);
%     end
% end

% for j = 1:n_terms
%     jdx = global2local(j,nsub);
%     for i = 1:n_terms
%         idx = global2local(i,nsub);
%         % product of A terms
%         % Atmp = tmp(i);
%         Atmp = sym(1);
%         for k = 1:total_order
%             ii = idx(k);
%             jj = jdx(k);
%             Atmp = Atmp*A(ii,jj);
%         end
%         dp = dp + Atmp*tmp(i);
%     end
% end
% y = A*x + b;
% dp = subs(dp,sym2cell(x),sym2cell(y));
end




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

function iG = local2global(iSub,nSub)
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