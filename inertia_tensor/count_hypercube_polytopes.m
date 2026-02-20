function counts = count_hypercube_polytopes(N)
    counts = zeros(N+1,1);
    for i = 0:N
        counts(i+1) = 2^(N-i)*nchoosek(N,i);
    end
end