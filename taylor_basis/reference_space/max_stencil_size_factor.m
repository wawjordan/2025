function [rat,max_sz,min_sz] = max_stencil_size_factor(n_dim,degree)
min_sz = ceil(1.5*nchoosek( n_dim + degree, degree ));
szs    = crystal_ball_sequence(n_dim,0:degree+2);
loc = find(szs-min_sz>0,1,"first");
max_sz = szs(loc);

rat = ceil(max_sz/min_sz);

end

% f = @(n_dim,degree) arrayfun(@(degree) crystal_ball_sequence(n_dim,degree)./nchoosek( n_dim + degree, degree )./(2.^(n_dim)), degree );