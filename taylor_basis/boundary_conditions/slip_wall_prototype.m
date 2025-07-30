%% prototype script for slip-wall BC for variational rec. (06/22/2025)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = '2025';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 3;
neq = 1;
syms u [d 1]
syms n [d 1]
syms t [d 1]
syms s [d 1]
v = [n.';t.';s.'];
v = v(1:d,1:d);
% M1 = kron(v(1:d,1:d),eye(neq));
% v2 = [((-1).^(1:neq).');repmat(((-1).^(0:neq-1).'),d-1,1)];
% M2 = diag(v2) * M1;
% A1 = M2\M1;
Dp = eye(d) - 2*diag(unit_vec(d,1));

A = v.' * Dp * v;
A2 = v.' * v - 2*n*(n.');
A3 = eye(d) - 2*n*(n.');
simplify( A == A2 )

M = @(p) (-1)^p * A3 - eye(d);
% A2 = @(p) v.'*diag( ones(d,1) - (1 +(-1)^p) * unit_vec(d,1) )*v - eye(d);
% A2_alt = @(p) -(-1)^p * n*(n.') + s*(s.') + t*(t.') - eye(d);
% 
% A2_alt2 = @(p) v.'*v -(1 +(-1)^p)*n*(n.') - eye(d);
% simplify( A2(2) == A2_alt(2) )


function e = unit_vec(d,n)
I = eye(d);
e = I(:,n);
end

% vec = @(x) x(:);
% simplify( inv(v)*det(v) == [ cross( v(2,:), v(3,:) ).', ...
%                              cross( v(3,:), v(1,:) ).', ...
%                              cross( v(1,:), v(2,:) ).' ] )
% thus inv(v) == v.'

% [ cross(t,s), cross(s,n), cross(n,t) ]

% for i=1:3
%     I = circshift(eye(d),i-1,1);
%     A_1 = subs( A1, num2cell( vec(symvar(A1)) ), num2cell( I(:) ) );
%     A_2 = subs( A2, num2cell( vec(symvar(A2)) ), num2cell( I(:) ) );
%     all(A_1==A_2,"all")
% end
% A2 = @(p) v.'*diag( ones(d,1) - (1 +(-1)^p) * unit_vec(d,1) )*v - eye(d);
% A2_alt = @(p) -(-1)^p * n*(n.') + s*(s.') + t*(t.') - eye(d);
% 
% A2_alt2 = @(p) v.'*v -(1 +(-1)^p)*n*(n.') - eye(d);
% simplify( A2(2) == A2_alt(2) )