function [v1,v2] = multichoose(n,k)

% def 1
v1 = factorial(n)/prod( factorial(k) );

% def 2
v2 = 1;
m = numel(k);
for i = 1:m
    v2 = v2 * nchoosek( sum(k(1:i)), k(i) );
end

end