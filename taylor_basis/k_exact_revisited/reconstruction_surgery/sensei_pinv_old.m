function LHS_plus = sensei_pinv_old(LHS)
[U,S,V] = svd(LHS,'vector');
% LHS = U*diag(S)*V'
% integer,                          intent(in)  :: LHS_m, LHS_n
% real(dp), dimension(LHS_m,LHS_n), intent(in)  :: LHS
% real(dp), dimension(LHS_n,LHS_m), intent(out) :: LHS_plus
% Compute Inverse of Singular Values
LHS_m = size(LHS,1);
LHS_n = size(LHS,2);
LHS_plus = zeros(LHS_n,LHS_m);
min_LHS_extents = min(LHS_m,LHS_n);

Sinv = zeros(min_LHS_extents,1);

% Determine Truncation Tolerance
% machine_precision = (10)**(-precision(one))
machine_precision = 1e-15;
abs_tolerance     = sqrt(machine_precision);
rel_tolerance     = S(1)*abs_tolerance;

for i = 1:min_LHS_extents
    if ( S(i) <= rel_tolerance )
        % Truncate Singular Value
        Sinv(i) = 0;
    else
        Sinv(i) = 1/S(i);
    end
end

% Compute Pseudo-Inverse
for i = 1:LHS_n
    V_Sinv = zeros(LHS_m,1);
    for j = 1:LHS_n
        % V_Sinv(j) = VT(j,i)*Sinv(j)
        V_Sinv(j) = V(i,j)*Sinv(j);
    end
    for j = 1:LHS_m
        LHS_plus(i,j) = dot( V_Sinv, U(j,:) );
    end
end

end