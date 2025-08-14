function exp = get_exponents_2D(degree)
% N = nchoosek(2+degree,degree);
% exp = zeros(2,N);
cnt = 0;
for j = 0:degree
    for i = j:degree
        cnt = cnt + 1;
        exp(1,cnt) = j;
        exp(2,cnt) = i-j;
    end
end

end