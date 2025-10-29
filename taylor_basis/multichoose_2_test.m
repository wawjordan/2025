function multichoose_2_test(n_dim,max_N,n_iter)
if ( max_N > 10 )
    error('that is way too big')
end
alpha = zeros(n_dim,1);
beta = zeros(n_dim,1);
for i = 1:n_iter
    for j = 1:n_dim
        alpha(j) = randi(max_N);
        beta(j)  = randi(alpha(j));
    end
    val = multichoose_2(alpha,beta);
    if (val(1) ~= val(2))
        fprintf('%s\n',sprintf(sprintf('(%d%%d )',n_dim),alpha));
        fprintf('%s\n',sprintf(sprintf('(%d%%d )',n_dim),beta));
        fprintf('val(1)=%d\n',val(1));
        fprintf('val(2)=%d\n',val(2));
        error('values do not match')
    end
end

end

