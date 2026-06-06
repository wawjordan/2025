function [x,dx,t,xb,dxb,tb] = combine_curves(N_i,x_i,d_i)
% given a list of abscissae (ti) and spacings (di) compute a C1 stretching function
N_i = N_i(:); % node count at xi (Ni(1)=1)
x_i = x_i(:); % location in physical space at Ni
d_i = d_i(:); % spacing at control point
s_i = get_slopes(N_i,x_i,d_i);

refine = true;
n_sub = numel(N_i) - 1;

x   = zeros(N_i(end),1);
dx  = zeros(N_i(end),1);
t   = zeros(N_i(end),1);
xb  = zeros(n_sub+1,1);
dxb = zeros(n_sub+1,1);
tb  = zeros(n_sub+1,1);
tb(1) = 0;
xb(1) = 0;

for i = 1:n_sub
    L  = x_i(i+1)-x_i(i);
    % LN = ( x_i(i+1)-x_i(i) )/( N_i(i+1)-N_i(i) );
    LN = ( x_i(i+1)-x_i(i) );
    t_ = linspace(0,1,N_i(i+1)-N_i(i)+1);
    [x_,dx_,~] = vinokur_two_sided(t_,s_i(i,1),s_i(i,2),refine);
    t_  = (x_i(i) + L*t_)/x_i(end);
    x_  = (x_i(i) + L*x_);
    x( N_i(i):N_i(i+1)) =  x_;
    dx(N_i(i):N_i(i+1)) = dx_*LN;
    t( N_i(i):N_i(i+1)) =  t_;

    xb(i+1)  = x_(end);
    dxb(i+1) = dx_(end);
    tb(i+1)  = t_(end);
end
end