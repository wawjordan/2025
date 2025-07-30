%% Isentropic Euler Vortex Solution (07/30/2025)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
inputs     = sym('inputs',[4,1]);
r_i   = inputs(1);
rho_i = inputs(2);
p_i   = inputs(3);
M_i   = inputs(4);

ref     = sym('ref',[3,1]);
l_ref   = ref(1);
rho_ref = ref(2);
a_ref   = ref(3);

syms gamma x y z real
xg   = 1/gamma;
gm1  = gamma - 1;
xgm1 = 1/gm1;
speed_of_sound = @(p_,rho_) sqrt( gamma.*p_./rho_);

% non-dimensionalization
r_i   = r_i / l_ref;
rho_i = rho_i / rho_ref;
p_i   = p_i / (rho_ref * a_ref.^2);
a_i   = speed_of_sound(p_i,rho_i);
u_i   = M_i * rho_i.^(gm1/2);

%
r    = sqrt( x.^2 + y.^2);
rho(x,y,z)  = rho_i .* ( 1 + (gm1/2)*M_i.^2 .* (1 - (r_i./r).^2 ) ).^xgm1;
p(x,y,z)    = xg .* rho.^gamma;
uvel(x,y,z) =  y .* a_i .* u_i .* r_i ./ r.^2;
vvel(x,y,z) = -x .* a_i .* u_i .* r_i ./ r.^2;
wvel(x,y,z) = z;
u = {uvel;vvel;wvel};

matlabFunction(rho ,"File","dens_svf","Vars",{x,y,z,gamma,inputs,ref},'Optimize',false);
matlabFunction(u{1},"File","uvel_svf","Vars",{x,y,z,gamma,inputs,ref});
matlabFunction(u{2},"File","vvel_svf","Vars",{x,y,z,gamma,inputs,ref});
matlabFunction(u{3},"File","wvel_svf","Vars",{x,y,z,gamma,inputs,ref});
matlabFunction(p   ,"File","pres_svf","Vars",{x,y,z,gamma,inputs,ref});
