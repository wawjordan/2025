function P = generate_poly_fun(dim,total_degree,varargin)
[dim,total_degree,n_terms] = pre_parse_inputs(dim,total_degree);
[exp,idx,~] = get_exponents( dim, total_degree );
default_rand = false;
default_coefs = ones(1,n_terms);
default_shift = zeros(1,dim);
default_scale = ones(1,dim);
p = inputParser;
validScale           = @(x) all(isnumeric(x)) && numel(x)==dim && all(x > 0);
validShift           = @(x) all(isnumeric(x)) && numel(x)==dim;
validCoefs           = @(x) all(isnumeric(x)) && numel(x)==n_terms;
addOptional(p,'coefs',default_coefs,validCoefs);
addOptional(p,'rand_coefs',default_rand,@(x)islogical(x));
addOptional(p,'scale',default_scale,validScale);
addOptional(p,'shift',default_shift,validShift);

parse(p,varargin{:});


if p.Results.rand_coefs
    [c_num,c_den] = rat( rand(1,n_terms) );
    for i = 2:n_terms
        c_den(i) = c_den(i)*prod(factorial(exp(:,i)));
    end
    coefs = c_num./c_den;
else
    coefs = p.Results.coefs(:).';
end
shift = p.Results.shift(:).';
scale = p.Results.scale(:).';


P = struct();
P.dim          = dim;
P.total_degree = total_degree;
P.n_terms      = n_terms;
P.idx          = idx;
P.exponents    = exp;
P.coefs        = coefs;
P.shift        = shift;
P.scale        = scale;


syms x [1 dim]
syms xc [1 dim]
syms L [1 dim]
syms a [1 n_terms]

x_ = (x-xc)./L;

fun_expr1 = a(1);
for i = 2:n_terms
    tmp_var1 = a(i);
    for d = 1:dim
        if (exp(d,i) == 0)
            continue
        end
        tmp_var1 = tmp_var1 .* x(d).^exp(d,i);
    end
    fun_expr1 = fun_expr1 + tmp_var1;
end
% fun_expr = recursive_horner(fun_expr,x);
% recover the grevlex/grelex ordering
[~,ind] = sort(coeffs(fun_expr1,x));

% add in the shift and scale
fun_expr2 = subs(fun_expr1,x,x_);
out_coefs = coeffs(fun_expr2,x);
out_coefs = out_coefs(ind);
% matlabFunction(vpa(out_coefs(ind)),"Vars",{[a(:)],[xc(:)],[L(:)]})

tmp_expr = subs(fun_expr2,a,sym(coefs));
tmp_expr = subs(tmp_expr,xc,sym(shift));
fun(x) = subs(tmp_expr,L,sym(scale));
P.test_fun = matlabFunction(vpa(fun));

% now for the derivatives ...
P.deriv = cell(1,total_degree);

for n = 1:total_degree
    n_terms_start = nchoosek( dim + n - 1, n - 1);
    n_terms_end   = nchoosek( dim + n, n );
    n_terms_local = n_terms_end - n_terms_start;
    P.deriv{n} = cell(1,n_terms_local);
    for i = n_terms_start+1:n_terms_end
        tmp_var1 = fun_expr1;
        tmp_var2 = fun_expr2;
        for d = 1:dim
            if (exp(d,i) == 0)
                continue
            end
            for j = exp(d,i):-1:1
                tmp_var1 = diff(tmp_var1,x(d));
                tmp_var2 = diff(tmp_var2,x(d));
            end
        end
        % tmp_var = recursive_horner(tmp_var,x);
        % recover the grevlex/grelex ordering
        [~,ind] = sort(coeffs(tmp_var1,x));

        % add in the shift and scale
        out_coefs = coeffs(tmp_var2,x);
        out_coefs = out_coefs(ind);
        
        tmp_expr = subs(tmp_var2,a,sym(coefs));
        tmp_expr = subs(tmp_expr,xc,sym(shift));
        fun(x) = subs(tmp_expr,L,sym(scale));
        P.deriv{n}{i-n_terms_start} = matlabFunction(vpa(fun));
    end
end
P.d_eval = @(order,varargin) eval_derivative(P,order,varargin{:});
end

function val = eval_derivative(P,order,varargin)
  if ( numel(order) ~= P.dim )
      error('incorrect number of dimensions for derivative order')
  end
  if ( numel(varargin) ~= P.dim)
      error('incorrect number of space coordinates')
  end
  if ( any(order<0) )
      error('incorrect number of space coordinates')
  end
  if sum(order)==0 % just evaluate the function
      val = P.test_fun(varargin{:});
  else
      cnt = 1;
      for j = 1:P.total_degree
        for i = 1:P.idx(j+1)-P.idx(j)
            cnt = cnt + 1;
          if all(order(:)==P.exponents(:,cnt))
              % val = P.deriv{j}{i}(varargin{:});
              val = P.deriv{j}{i}(varargin{:}) + varargin{1}*0;
              return
          end
        end
      end
      % if it makes it here, then the result should be zero
      val = 0;
  end
end

function [dim_out,total_degree_out,n_terms_out] = pre_parse_inputs(dim,total_degree)
p = inputParser;
validScalarNonNegNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
addRequired(p,'dim',validScalarNonNegNum);
addRequired(p,'total_degree',validScalarNonNegNum);
parse(p,dim,total_degree);
dim_out = p.Results.dim;
total_degree_out = p.Results.total_degree;
n_terms_out = get_n_terms(dim_out,total_degree_out);
end

function expr_out = recursive_horner(expr_in,vars)
expr_out = expr_in;
for d = 1:numel(vars)
    expr_out = horner(expr_out,vars(d));
end
end

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