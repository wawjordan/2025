%% SVF test case with parameterized grid (07/24/2026)
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

clc; clear; close all;

folder = 'C:\Users\Will\Downloads\SVF_';
jobfmt = '_svf%0.4dx%0.4d';

% Nfine = [2049,1025];
% N_levels = 8;

% [N_theta,N_r,AR_out,r,radius] = find_N_cells(AR_target,min_N,max_N,N_refine,angle,radius1,radius2)
n_ref = 0;
[N_theta,N_r,AR_out,r,radius] = find_N_cells(1,33,257,n_ref,pi/2,2,3);
Ntheta_2 = (N_theta-1)/2^n_ref + 1;
Nr_2 = (N_r-1)/2^n_ref + 1;
GRID = svf_grid(N_theta,radius);
hold on
% plot(GRID.x(1:2^n_ref:end,1:2^n_ref:end),  GRID.y(1:2^n_ref:end,1:2^n_ref:end),  'k');
% plot(GRID.x(1:2^n_ref:end,1:2^n_ref:end).',GRID.y(1:2^n_ref:end,1:2^n_ref:end).','k');
% axis equal

AR0 = AR_out;
ARN = AR_out;
[radius2,~,AR_out2] = get_stretched_radius_distribution(AR0,ARN,N_r,N_theta,pi/2,2,3);
GRID2 = svf_grid(N_theta,radius2);
plot(GRID2.x(1:2^n_ref:end,1:2^n_ref:end),  GRID2.y(1:2^n_ref:end,1:2^n_ref:end),  'r');
plot(GRID2.x(1:2^n_ref:end,1:2^n_ref:end).',GRID2.y(1:2^n_ref:end,1:2^n_ref:end).','r');
axis equal
figure
plot(1./AR_out2)
hold off
% r = r_from_n_theta(AR,n_theta,angle);
% hold on
% for N = 16:32
%     fplot(@(r) arrayfun(@(r) aspect_ratio_from_d_theta(r,N,pi/2),r),[1,3],'g')
%     fplot(@(r) aspect_ratio_from_d_theta2(r,N,pi/2),[1,3],'r')
%     fplot(@(r) aspect_ratio_from_d_theta3(r,N,pi/2),[1,3],'b')
%     delta_theta = (pi/2)/(N-1);
%     plot((2+delta_theta)/(2-delta_theta),1,'o')
% end

% inner radius 2
% outer radius 3
% angle: pi/2
% arc length = r*(pi/2)
% inner side cell length: pi*(1/(Ni-1))
% outer side cell length: (3/2)*pi*(1/(Ni-1))
% How to distribute along j-direction to keep ~unit aspect ratio?



% for j = 1:N_levels
%     s = 2^(j-1);
%     GRID2 = grid_subset(GRID,{1:s:Nfine(1),1:s:Nfine(2),1});
%     foldername = [folder,sprintf(jobfmt,GRID2.imax,GRID2.jmax)];
%     filename = [foldername,'\svf'];
%     status = mkdir(foldername);
%     svf_bc_out(GRID2.imax,GRID2.jmax,[filename,'.bc'])
% end
% function [N,r,AR_out] = get_radius_distribution(AR_target,n_refine,n_theta,angle,radius1,radius2)
% delta_theta = angle/(n_theta-1);
% r = (2 + AR_target*delta_theta)/(2 - AR_target*delta_theta);
% n1 = log( r*((radius2-radius1)/radius1 + 1) )/log(r);
% res = (n1-1)/2^n_refine;
% N = floor(res)*2^n_refine + 1;
% [~,dx0,~] = my_geomspace(N,radius1,xmax=radius2,r=r);
% AR_out = ( dx0/(2*radius1+dx0) ) * (2/delta_theta);
% end

function [r,r_fun,AR_out] = get_stretched_radius_distribution(AR_0,AR_N,n_r,n_theta,angle,radius1,radius2)
delta_theta = angle./(n_theta-1);
alpha0 = (2 + AR_0.*delta_theta)./(2 - AR_0.*delta_theta);
alphaN = (2 + AR_N.*delta_theta)./(2 - AR_N.*delta_theta);
delta_r0 = (alpha0-1)*radius1;
delta_rN = (1-1/alphaN)*radius2;
delta_r  = radius2-radius1;
fun = vinokur_two_sided_spacing_fcn( n_r, delta_rN/delta_r, delta_r0/delta_r,false);
r_fun = @(t) radius1 + delta_r*fun(t);
r = r_fun(linspace(0,1,n_r));
alpha = r(2:end)./r(1:end-1);
AR_out = ( (alpha-1)./(alpha+1) ) .* (2./delta_theta);
end

function [N,r,AR_out] = get_radius_distribution(AR_target,n_refine,n_theta,angle,radius1,radius2)
delta_theta = angle./(n_theta-1);
r = (2 + AR_target.*delta_theta)./(2 - AR_target.*delta_theta);
n1 = log(radius2./radius1)./log(r);
res = (n1-1)./2.^n_refine;
N = ceil(res).*2.^n_refine + 1;
r = (radius2./radius1).^(1./(N-1));
AR_out = ( (r-1)./(r+1) ) .* (2./delta_theta);
end

function [N_theta,N_r,AR_out,r,radius] = find_N_cells(AR_target,min_N,max_N,N_refine,angle,radius1,radius2)
sz_check = @(x) mod(x-1,2.^N_refine)==0 & (x-1)./(2.^N_refine)>=(min_N-1) & x<=max_N;
N_in = min_N:max_N;
N_in = N_in( sz_check(N_in) );
[N1,r1,AR] = get_radius_distribution(AR_target,N_refine,N_in,angle,radius1,radius2);
tol = 1.0e-14;
mask = sz_check(N1) & abs(radius2./radius1 - r1.^(N1-1))<tol;
N_in = N_in(mask);
N1   = N1(mask);
r1   = r1(mask);
AR   = AR(mask);
[~,idx] = min( abs(AR-AR_target) );

N_theta = N_in(idx);
N_r     = N1(idx);
AR_out  = AR(idx);
r       = r1(idx);
radius = r.^(0:N_r-1).*radius1;
end


function AR = aspect_ratio_from_d_theta(r,n_theta,angle)
delta_theta = angle/(n_theta-1);
x = zeros(4,1);
y = zeros(4,1);

x(1) = cos( delta_theta/2);
x(2) = cos(-delta_theta/2);
x(3) = r*cos(-delta_theta/2);
x(4) = r*cos( delta_theta/2);

y(1) = sin( delta_theta/2);
y(2) = sin(-delta_theta/2);
y(3) = r*sin(-delta_theta/2);
y(4) = r*sin( delta_theta/2);
% AR = aspect_ratio_2(x,y);
AR = aspect_ratio_1(x,y);
end

function AR = aspect_ratio_from_d_theta2(r,n_theta,angle)
delta_theta = angle./(n_theta-1);
% AR = ((r-1)./(r+1))./(delta_theta/2);
AR = ( 2./(r+1) ).*((r-1)./delta_theta);
AR = min(AR,1./AR);
end

function AR = aspect_ratio_from_d_theta3(r,n_theta,angle)
delta_theta = angle./(n_theta-1);
rc = (4/3)*(sin(delta_theta/2)./delta_theta).*(r.^3-1)./(r.^2-1);
AR = ((r-1)./rc)./(delta_theta);
AR = min(AR,1./AR);
end

function val = aspect_ratio_1(x,y)
m1  = sqrt( ( ( x(2) - x(1) ) + (x(3) - x(4) ) )^2 ...
          + ( ( y(2) - y(1) ) + (y(3) - y(4) ) )^2 );
m2  = sqrt( ( ( x(3) - x(2) ) + (x(4) - x(1) ) )^2 ...
          + ( ( y(3) - y(2) ) + (y(4) - y(1) ) )^2 );
val = m1/m2;
% val = max( val, 1/val );
val = min( val, 1/val );
end

function val = aspect_ratio_2(x,y)
hmax = max( sqrt( (x(3)-x(1))^2 + (y(3)-y(1))^2 ), ...
            sqrt( (x(4)-x(2))^2 + (y(4)-y(2))^2 ) );
hmag = 0;
minAx2 = realmax;
for k = 1:4
    kp1  = mod(k  ,4)+1;
    kp3  = mod(k+2,4)+1;
    minAx2 = min( minAx2, ...
                ( x(k)   - x(kp3) )*( y(kp1) - y(kp3)) ...
              - ( y(k)   - y(kp3) )*( x(kp1) - x(kp3) ) );
    h2   = (x(kp1)-x(k))^2 + (y(kp1)-y(k))^2;
    hmag = hmag + h2;
    hmax = max( hmax, sqrt(h2) );
end
% val = 0.25*sqrt(2*hmag)*hmax/minAx2;
val = 1/( 0.25*sqrt(2*hmag)*hmax/minAx2 );
end

function GRID = svf_grid(N_theta,radius)

GRID.dim  = 2;
N_z = 1;
GRID.imax = N_theta;
GRID.jmax = numel(radius);
GRID.kmax = N_z;
t_range = [pi/2, 0];
z_range = [0, 0];

t = linspace(t_range(1),t_range(2),N_theta);
z = linspace(z_range(1),z_range(2),N_z);

[T,R,GRID.z] = ndgrid(t,radius,z);

GRID.x = R.*cos(T);
GRID.y = R.*sin(T);

end

% function svf_bc_out(N1,N2,filename)
% fid = fopen(filename,'w');
% intfmt = '        %-5d';
% % fltfmt = ' %-# 23.16E';
% 
% fprintf(fid,'%d\n',1);
% fprintf(fid,'%d\n',1);
% fprintf(fid,'%d     %d\n',N1,N2);
% fprintf(fid,'A\n');
% fprintf(fid,'%d\n',4);
% fprintf(fid,['bc1   ',repmat(intfmt,[1,5]),'\n'], 1,N1, 1, 1,-200);
% fprintf(fid,['bc2   ',repmat(intfmt,[1,5]),'\n'],N1,N1, 1,N2, 601);
% fprintf(fid,['bc3   ',repmat(intfmt,[1,5]),'\n'], 1,N1,N2,N2,-200);
% fprintf(fid,['bc4   ',repmat(intfmt,[1,5]),'\n'], 1, 1, 1,N2,-200);
% 
% fclose(fid);
% end
% 
% function svf_bc_out_mms(N1,N2,filename)
% fid = fopen(filename,'w');
% intfmt = '        %-5d';
% % fltfmt = ' %-# 23.16E';
% 
% fprintf(fid,'%d\n',1);
% fprintf(fid,'%d\n',1);
% fprintf(fid,'%d     %d\n',N1,N2);
% fprintf(fid,'A\n');
% fprintf(fid,'%d\n',4);
% fprintf(fid,['bc1   ',repmat(intfmt,[1,5]),'\n'], 1,N1, 1, 1,-200);
% fprintf(fid,['bc2   ',repmat(intfmt,[1,5]),'\n'],N1,N1, 1,N2,-200);
% fprintf(fid,['bc3   ',repmat(intfmt,[1,5]),'\n'], 1,N1,N2,N2,-200);
% fprintf(fid,['bc4   ',repmat(intfmt,[1,5]),'\n'], 1, 1, 1,N2,-200);
% 
% fclose(fid);
% end