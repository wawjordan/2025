function coef_struct = get_linear_reconstruction_of_function_over_stencil_stencil(TS,QS,fun)
n_nbors = numel(TS);

coef_struct = struct();
coef_struct.coefs = cell(n_nbors,1);
for n1 = 1:n_nbors
    % get cell averages
    avgs = get_cell_avg_of_function(QS(n1).Q,fun);
    coef_struct(n1).avgs = avgs;

    % now go through and compute the linear reconstructions
    nn_nbors = numel(TS(n1).T);

    M1 = [];
    R1 = [];
    for n2 = 2:nn_nbors
        M1 = [ M1; TS(n1).T(1).build_linear_reconstruction_matrix( TS(n1).T(n2) ) ];
        R1 = [ R1; TS(n1).T(1).build_linear_reconstruction_rhs( TS(n1).T(n2), avgs(1), avgs(n2) )];
    end

    coefs = [avgs(1);M1\R1];

    coef_struct(n1).coefs = coefs;

end