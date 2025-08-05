%% combinatorics puzzle (08/05/2025)
clc; clear; close all;

dim    = 1;
degree = 3;
n_terms = get_n_terms_( dim, degree );
order = get_exponents_sorted_( dim, degree );

% num_fmt = '%4d';
num_fmt = '%5d/%-5d';

fmt_fun = @(n,fmt) [repmat([fmt,' '],1,n-1), fmt, '\n'];

fmt = fmt_fun(4,num_fmt);

for i = 1:n_terms
    den = get_factorial_scaling_1( order(:,i) );
    num1 = get_factorial_scaling_2( order(:,i) );
    [num2,num3] = multichoose( sum(order(:,i)), order(:,i) );
    tmp = 1/prod( factorial(order(:,i)));
    [N,D] = rat([num1/den,num2/den,num3/den,tmp]);
    tmp_vec = [N;D];
    fprintf(fmt,tmp_vec(:));
    % fprintf(fmt,[den1,den2]);
end


% num3 = nchoosek(sum(order(1:n_dim)),order(1)) * nchoosek(sum(order(2:n_dim)),order(2));

function den = get_factorial_scaling_1(order)
den = factorial( sum(order) );
end

function num = get_factorial_scaling_2(order)
n_dim = numel(order);
num = 1;
for d = 1:n_dim
    num = num * nchoosek( sum( order(d:n_dim) ), order(d) );
end
end


function n_terms = get_n_terms_( dimen, total_degree )
if (total_degree < 0)
    n_terms = 0;
else
    n_terms = nchoosek( dimen + total_degree, total_degree );
end
end

function exp = get_exponents_sorted_( dimen, total_degree )
            n_terms = get_n_terms_(dimen,total_degree);
            exp     = zeros(dimen,n_terms);
            p_idx     = zeros(n_terms,1);
            full_degree = (total_degree+1)^dimen;
            cnt = 0;
            nSub = (total_degree+1)*ones(dimen,1);
            n10s = 10^floor(log10(full_degree));
            for n = 0:full_degree
                tmp_exp = global_to_local_(n+1,nSub) - 1;
                current_degree = sum(tmp_exp);
                if current_degree <= total_degree
                    cnt = cnt + 1;
                    exp(:,cnt) = tmp_exp;
                    % (order,reverse lexicographic in each dim) sort
                    p_idx(cnt) = current_degree*n10s*10^(dimen-1);
                    for j = 1:dimen
                        p_idx(cnt) = p_idx(cnt) + exp(j,cnt)*10^(j-1);
                    end

                end
            end
            [~,idx] = sort(p_idx);
            exp = exp(:,idx);
end


function iSub = global_to_local_(iG,nSub)
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