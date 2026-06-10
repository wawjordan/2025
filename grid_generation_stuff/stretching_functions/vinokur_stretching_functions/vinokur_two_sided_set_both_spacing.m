function [s0,s1] = vinokur_two_sided_set_both_spacing(N,d0,d1,refine)

d = [d0,d1];
dt = 1/(N-1);
S =d/dt;
lb = [1e-6,1e-6];
ub = [1e6,1e6];
options = optimoptions("fmincon",...
                       "Algorithm","interior-point", ...
                       "HonorBounds",true, ...
                       "FunctionTolerance",1e-16,...
                       "OptimalityTolerance",1e-16,...
                       "StepTolerance",1e-30,...
                       "Display","none");
s = fmincon(@(s)obj_fun(N,s,d,refine),S,[],[],[],[],lb,ub,[],options);
% options = optimoptions("fminunc","FunctionTolerance",1e-16,"OptimalityTolerance",1e-16,"Display","iter");
% s = fminunc(@(s)obj_fun(N,s,d),S,options);
s0 = 1/s(1);
s1 = 1/s(2);
end

function f = obj_fun(N,s,d0,refine)
    d = vinokur_two_sided_spacing_start_end(N,s(1),s(2),refine);
    f = sum((d(:)-d0(:)).^2);
end

function d = vinokur_two_sided_spacing_start_end(N,s0,s1,refine)
dt = 1/(N-1);
[d,~,~] = vinokur_two_sided([dt,1-dt],s0,s1,refine);
d = d(:);
d(2) = 1-d(2);
end