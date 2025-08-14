function [exp,val] = eval_basis_2D_at_quad_pts(degree,x0,h,xq,wq)
N = nchoosek(2+degree,degree);
nquad = numel(wq);

exp = zeros(2,N);
val = ones(nquad,N);
dq = zeros(nquad,2);
mn = zeros(N,1);

xh = 1./h(:);

for n = 1:nquad
    dq(n,1) = (xq(1,n) - x0(1)).*xh(1);
    dq(n,2) = (xq(2,n) - x0(2)).*xh(2);
end

cnt = 0;
for j = 0:degree
    for i = 0:degree-j
        cnt = cnt + 1;
        exp(1,cnt) = i;
        exp(2,cnt) = j;
        for ii = 1:i
            val(:,cnt) = val(:,cnt).*dq(:,1);
        end
        for jj = 1:j
            val(:,cnt) = val(:,cnt).*dq(:,2);
        end

        for n = 1:nquad
            mn(cnt) = mn(cnt) + val(n,cnt).*wq(n);
        end
    end
end

vol = sum( wq );

mn = mn / vol;

for cnt = 1:N
    val(:,cnt) = val(:,cnt) - mn(cnt);
end

end