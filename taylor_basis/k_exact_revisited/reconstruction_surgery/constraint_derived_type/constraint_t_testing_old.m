%% constraint_t_old testing (06/04/2024)
clc; clear; close all;

% poly_degree = 2;
% n_var       = 2;
% n_dim       = 3;
% 
% p = polynomial_mod_t_old(n_var, n_dim, poly_degree);
% p.coeff = rand(p.n_var,p.n_terms);
% p.x_ref = rand(p.n_dim,1);
% constraint_type = 1000;
% constraint_eval = 1;
% constraint_inputs = constraint_inputs_t_old();
% constraint_inputs.VAR_IDX = 1;
% constraint_inputs.x_eval(:) = 0;
% constraint_inputs.x_ref(:)  = 0;
% constraint_inputs.normal(:) = [1;0;0];
% constraint_inputs.coupling_coeff(:) = 0;
% 
% 
% con1 = constraint_t_old( poly_degree, n_var, n_dim, constraint_type, constraint_eval, constraint_inputs );
% con1 = con1.generate_constraint(p);
% con1.constraint_row

poly_degree = 2;
n_var       = 5;
n_dim       = 3;
p = polynomial_mod_t_old(n_var, n_dim, poly_degree);
p.coeff = rand(p.n_var,p.n_terms);
p.x_ref = rand(p.n_dim,1);
normal = rand(p.n_dim,1); normal = normal/norm(normal,2);

x_eval = rand(p.n_dim,1);

constraint_type = 3000;
constraint_eval = 0;

constraint_inputs = constraint_inputs_t_old();
constraint_inputs.x_eval(:)  = x_eval;
constraint_inputs.x_ref(:)   = p.x_ref;
constraint_inputs.normal(:)  = normal;
constraint_inputs.VAR_IDX = [2;3;4];
constraint_inputs.coupling_coeff = normal;

n_constraints = 1;
bc_constraints = constraint_container_t_old();
bc_constraints.n_constraints = n_constraints;
bc_constraints.constraint = repmat(constraint_t_old(),[n_constraints,1]);
for ci = 1:n_constraints
    bc_constraints.constraint(ci) = constraint_t_old( poly_degree, n_var, n_dim, constraint_type(ci), constraint_eval(ci), constraint_inputs(ci) );
    con1 = bc_constraints.constraint(ci).generate_constraint(p);
    con1.constraint_row
end

% con1 = constraint_t_old( poly_degree, n_var, n_dim, constraint_type, constraint_eval, constraint_inputs );
% con1 = con1.generate_constraint(p);
% con1.constraint_row