function val = multichoose_2(alpha,beta)

n1 = prod(factorial(alpha));
d1 = prod(factorial(beta));
d2 = prod(factorial(alpha-beta));
val(1) = n1 / (d1*d2);

val(2) = 1;
d = numel(alpha);
for i = 1:d
      val(2) = val(2) * nchoosek( alpha(i), beta(i) );
end

end