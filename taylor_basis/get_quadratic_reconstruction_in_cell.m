function coefs = get_quadratic_reconstruction_in_cell(Ts,Cs)

n_nbors = numel(Ts);

M1 = [];
R1 = [];
for n1 = 2:n_nbors
    M1 = [ M1; Ts(1).build_reconstruction_matrix( Ts(n1) ) ];
    R1 = [ R1; Ts(1).build_reconstruction_rhs(    Ts(n1), Cs(1).coefs, Cs(n1).coefs )];
end

coefs = [Cs(1).coefs;M1\R1];

end