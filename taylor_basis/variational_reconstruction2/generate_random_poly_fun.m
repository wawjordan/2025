function [test_fun,deriv] = generate_random_poly_fun(dim,total_degree)

exp = get_exponents_sorted( dim, total_degree );
n_terms = size(exp,2);
[c_num,c_den] = rat( rand(1,n_terms) );
for i = 2:n_terms
    c_den(i) = c_den(i)*prod(factorial(exp(:,i)));
end
% c_num = 0*c_num + 1;
% c_den = 0*c_den + 1;
syms x [1 dim]
syms a [1 n_terms]
fun_expr = a(1);
for i = 2:n_terms
    tmp_var = a(i);
    for d = 1:dim
        if (exp(d,i) == 0)
            continue
        end
        tmp_var = tmp_var .* x(d).^exp(d,i);
    end
    fun_expr = fun_expr + tmp_var;
end
% test_fun1 = matlabFunction(vpa(fun));
% for d = 1:dim
%     fun_expr = horner(fun_expr,x(d));
% end
fun_expr = recursive_horner(fun_expr,x);
fun(x) = subs(fun_expr,a,sym(c_num./c_den));
test_fun = matlabFunction(vpa(fun));

% now for the derivatives ...
deriv = cell(1,total_degree);

for n = 1:total_degree
    n_terms_start = nchoosek( dim + n - 1, n - 1);
    n_terms_end = nchoosek( dim + n, n );
    n_terms_local = n_terms_end - n_terms_start;
    deriv{n} = cell(1,n_terms_local);
    for i = n_terms_start+1:n_terms_end
        tmp_var = fun_expr;
        % exp(:,i)
        for d = 1:dim
            if (exp(d,i) == 0)
                continue
            end
            for j = exp(d,i):-1:1
                tmp_var = diff(tmp_var,x(d));
            end
        end
        tmp_var = recursive_horner(tmp_var,x);
        fun(x) = subs(tmp_var,a,sym(c_num./c_den));
        deriv{n}{i-n_terms_start} = matlabFunction(vpa(fun));
    end
end
end

function expr_out = recursive_horner(expr_in,vars)
expr_out = expr_in;
for d = 1:numel(vars)
    expr_out = horner(expr_out,vars(d));
end
end

function exp = get_exponents_sorted( dimen, total_degree )
n_terms = nchoosek( dimen + total_degree, total_degree );
exp     = zeros(dimen,n_terms);
p_idx     = zeros(n_terms,1);
full_degree = (total_degree+1)^dimen;
cnt = 0;
nSub = (total_degree+1)*ones(dimen,1);
n10s = 10^floor(log10(full_degree));
for n = 0:full_degree
    tmp_exp = global_to_local(n+1,nSub) - 1;
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