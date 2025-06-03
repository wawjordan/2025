function coef_struct = get_exact_linear_rec_of_function_over_stencil_stencil(TS,QS,fun,gradfun)
n_nbors = numel(TS);

coef_struct = struct();
coef_struct.coefs = cell(n_nbors,1);
for n1 = 1:n_nbors
    % get cell averages
    avgs = get_cell_avg_of_function(QS(n1).Q,fun);
    coef_struct(n1).avgs = avgs;
    % evaluate gradient at cell centroid
    xc = num2cell(TS(n1).T(1).x_ref);
    df = gradfun(xc{:});
    dx = TS(n1).T(1).h_ref;
    coefs = [avgs(1); df(:).*dx(:)];
    coef_struct(n1).coefs = coefs;
end

end