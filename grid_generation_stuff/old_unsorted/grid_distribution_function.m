function F = grid_distribution_function(y,d,N)
N_regions = length(N);
Fx = struct();
for i = 1:N_regions
    Fx(i).F = tanh_stretching_function(y(i),y(i+1),d(i),d(i+1),N(i));
end

F1 = @(x) arrayfun( @(i) Fx(i).F(x), discretize(bound(x,y(1),y(end)),y) );
F = @(x) arrayfun( @(x)F1(x), x );
end

function y = bound(x,bl,bu)
  % return bounded value clipped between bl and bu
  y=min(max(x,bl),bu);
end