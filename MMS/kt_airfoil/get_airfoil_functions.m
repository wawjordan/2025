function funs = get_airfoil_functions(epsilon,kappa,tau,alpha,scale)

vinf    = 75.0;
rhoinf  = 1.0;
pinf    = 100000.0;
gamma   = 1.4;

rho_ref = 1.0;
p_ref   = 100000.0;
a_ref   = sqrt(gamma*p_ref/rho_ref);

% nondimensionalize inputs
vinf   = vinf/a_ref;
rhoinf = rhoinf/rho_ref;
pinf   = pinf/(rho_ref*a_ref^2);

airfoil        = kt_airfoil(epsilon,kappa,tau);
airfoil.vinf   = vinf;
airfoil.rhoinf = rhoinf;
airfoil.pinf   = pinf;
airfoil        = airfoil.set_alpha(alpha);

funs = cell(5,1);

for v = 1:5
    funs{v} = @(x,y) arrayfun( @(x,y) airfoil.get_prim_variables(v,x,y,scale),x,y);
end

end