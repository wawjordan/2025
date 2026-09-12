function check_tensor_diff(T,f,xvec)
% function to check all entries in tensor derivative returned by tensor_diff
% inputs:
%  T     - multidimensional array of symbolic variables containing partial
%          derivatives of function f(x) of order n
%  f     - symbolic multivariate function
%  xvec  - vector of symbolic variables (xvec = [x(1),x(2),...])
nsub = size(T,ndims(T));
if nsub>1
    nsub = size(T);
else
    nsub = numel(T);
end
for i = 1:numel(T)
    idx = global2local(i,nsub);
    xs = sym2cell(xvec(idx));
    cidx = num2cell(idx);
    check = logical(simplify(T(cidx{:})==diff(f,xs{:})));
    if ~check
        fmt = ['[%f',repmat(',%f',1,numel(idx)-1),']'];
        s = sprintf('entry %s does not match',fmt);
        error(s,idx);
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