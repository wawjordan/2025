function T = tensor_diff(f,xvec,n)
% recursively calculates symbolic derivative tensor of multivariate function f
% inputs:
%  f     - symbolic multivariate function
%  xvec  - vector of symbolic variables (xvec = [x(1),x(2),...])
%  n     - integer specifying the derivative order
%      0 - the function itself; rank 0 tensor (scalar)
%      1 - gradient;            rank 1 tensor
%      2 - hessian;             rank 2 tensor
%     ...
% outputs:
%  T     - multidimensional array of symbolic variables containing partial
%          derivatives of order n
n_dim = numel(xvec);
if n == 0
    T = f;
    return
elseif n == 1
    T = sym( ones(n_dim,1) );
    for d = 1:n_dim
        T(d,1) = diff(f,xvec(d),1);
    end
    return
else
    T1 = tensor_diff(f,xvec,n-1);
    T  = diff(T1,xvec(1));
    for d = 2:n_dim
        T  = cat(n,T,diff(T1,xvec(d)));
    end
    return
end
end