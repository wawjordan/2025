function coefs = get_exact_quadratic_reconstruction_in_cell(Ts,Qs,fun,gradfun,hessfun)
avg = get_cell_avg_of_function(Qs(1),fun);

% evaluate gradient and hessian at cell centroid
xc = num2cell(Ts(1).x_ref);
dx = Ts(1).h_ref;


df1 = gradfun(xc{:});
% df1 = df1([2;1]);
df1 = df1([1;2]);
dx1 = dx([1;2]);

df2 = hessfun(xc{:});
% df2 = df2([4;2;1]);
df2 = df2([1;2;4]);
dx2 = [ dx(1)^2; dx(1)*dx(2); dx(2)^2; ];


% coefs = [avg; df1(:); df2(:) ];
coefs = [avg; df1(:).*dx1(:); df2(:).*dx2(:) ];


end