function si = get_slopes(N_i,x_i,d_i)
% given a list of abscissae (ti) and spacings (di) compute a C1 stretching function
N_i = N_i(:); % node count at xi (Ni(1)=1)
x_i = x_i(:); % location in physical space at Ni
d_i = d_i(:); % spacing at control point

n_sub = numel(N_i) - 1;

si = zeros(n_sub,2);
refine = true;

L = diff(x_i);
N = diff(N_i)+1;
dx = L./N;
d1 = d_i(1:end-1)./L;
d2 = d_i(2:end)./L;

% 1st subdomain
[si(1,1),si(1,2)] = vinokur_two_sided_set_both_spacing(N(1),d1(1),d2(1),refine);

if (n_sub>1)
    if (mod(n_sub,2)==0)

        % last subdomain
        [si(n_sub,1),si(n_sub,2)] = vinokur_two_sided_set_both_spacing(N(n_sub),d1(n_sub),d2(n_sub),refine);

        % work up from first to middle
        for i = 2:n_sub/2
            si(i,1) = si(i-1,2);
            s1 = si(i,1)/dx(i);
            [si(i,2),~] = vinokur_two_sided_set_one_spacing(N(i),d2(i),s1,refine);
        end

        % save the middle slope
        stmp1 = si(n_sub/2,2);

        % work back from last to middle
        for i = n_sub-1:-1:n_sub/2+1
            si(i,2) = si(i+1,1);
            s1 = si(i,2)/dx(i);
            [si(i,1),~] = vinokur_two_sided_set_one_spacing(N(i),d2(i),s1,refine);
        end

        % average the middle slope
        stmp2 = si(n_sub/2+1,1);
        savg = (stmp1+stmp2)/2;
        si(n_sub/2+1,1) = savg;
        si(n_sub/2,2)   = savg;
    else % odd
        % last subdomain
        [si(n_sub,1),si(n_sub,2)] = vinokur_two_sided_set_both_spacing(N(n_sub),d1(n_sub),d2(n_sub),refine);

        % work up until the middle subdomain
        for i = 2:(n_sub-1)/2
            si(i,1) = si(i-1,2);
            s1 = si(i,1)/dx(i);
            [si(i,2),~] = vinokur_two_sided_set_one_spacing(N(i),d2(i),s1,refine);
        end

        % work back until the middle subdomain
        for i = n_sub-1:-1:(n_sub+1)/2+1
            si(i,2) = si(i+1,1);
            s1 = si(i,2)/dx(i);
            [si(i,1),~] = vinokur_two_sided_set_one_spacing(N(i),d2(i),s1,refine);
        end

        % fill in the middle
        si((n_sub-1)/2+1,1) = si((n_sub-1)/2,2);
        si((n_sub-1)/2+1,2) = si((n_sub+1)/2+1,1);
    end
end
end