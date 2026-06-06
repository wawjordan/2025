% function [xi,dxi,ddxi] = blend_piecewise_linear(t,N)
function [xi,xi_f] = blend_piecewise_linear(n_sub,N_i,xi_i)
% n_sub : number of subdivisions
% Ni    : (n_sub)   number of points in a subdivision (includes endpoints)
% xi_i  : (n_sub+1) monotonic locations in [0,1] (include endpoints!)
% ol    : (n_sub-1) overlap at subdivision juncture for blending
% xi    : output

% total number of points (accounts for redundant endpoints)
Nt = sum(N_i(:)-1)+1;

xi  = zeros(Nt,1);
cnt = 1;
for j = 1:numel(Ni)
    dxi = (xi_i(j+1)-xi_i(j))/(Ni(j)-1);
    for i =2:Ni(j)
        cnt = cnt + 1;
        xi(cnt) = xi(cnt-1) + dxi;
    end
end

t_i = [0;cumsum(N_i-1)]+1;

xi_i
% p
t
% t1 locs
% sum j from 1 to i: (Nj-1)/(Nt-1) == 1

% t2 locs
% specified t locations





% dt = [0;( diff([ti(:);0])./Ni(:) )*sum(Ni(:))];
% breaks = diff(t)
% mask = (t)
% xi0 = 
% [fp,dfp,ddfp] = c2_blend_poly(t);
end