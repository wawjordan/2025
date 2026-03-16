function val = crystal_ball_sequence(n_dim,n)
val = arrayfun( @(n) A008288(n_dim,n-1), n );
end

function val = A008288(n,k)
% https://oeis.org/A008288
% val = 0;
% for d = 0:min(n,k)
%     val = val + nchoosek(k, d)*nchoosek(n+k-d, k);
% end
val = sum(arrayfun(@(d)nchoosek(k,d).*nchoosek(n+k-d,k),0:min(n,k)));
end