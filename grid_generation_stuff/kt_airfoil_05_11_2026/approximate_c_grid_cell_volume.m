function V = approximate_c_grid_cell_volume(x_jmin,y_jmin,delta,boundary_distance,jmax)

imax = numel(x_jmin);
% Inner boundary: x_jmin, y_jmin (ordered clockwise)
Gx1 = x_jmin;
Gy1 = y_jmin;

% outflow1
Gx2 = [x_jmin(end),x_jmin(end)];
Gy2 = [y_jmin(end),boundary_distance];

% outer boundary
[Gx3,Gy3] = body_pts(linspace(0,1,imax),boundary_distance,boundary_distance);

% outflow2
Gx4 = [x_jmin(1),x_jmin(1)];
Gy4 = [-boundary_distance,y_jmin(1)];

t = transfinite_interpolant_2D( G1x=Gx1,G1y=Gy1,G2x=Gx2,G2y=Gy2,...
                                G3x=Gx3,G3y=Gy3,G4x=Gx4,G4y=Gy4 );

xi_vec = linspace(-1,1,imax);
dx0 = delta/(2*boundary_distance);
% eta_vec = my_geomspace(jmax,-1,dx0=dx0,xmax=1);
eta_vec = linspace(-1,1,jmax);
[X,Y] = t.eval_grid(xi_vec,eta_vec);
V = 0;
end

function [x,y] = body_pts(t,r,l)
arc_length = pi*r + 2*l;
tc1 = l/arc_length;
tc2 = 1 - tc1;
x = zeros(size(t));
y = zeros(size(t));
mask1 = t<tc1;
mask2 = (t>=tc1)&(t<=tc2);
mask3 = t>tc2;

x(mask1) = l - l*t(mask1)/tc1;
y(mask1) = -r;

x(mask2) = r*cos(pi*(0.5+(tc2-t(mask2))/(tc2-tc1)) );
y(mask2) = r*sin(pi*(0.5+(tc2-t(mask2))/(tc2-tc1)) );

x(mask3) = l*(t(mask3)-tc2)/(1-tc2);
y(mask3) =  r;
end