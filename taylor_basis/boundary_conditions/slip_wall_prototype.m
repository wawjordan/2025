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
M1 = kron(v(1:d,1:d),eye(neq));
% M1 = [ kron(n.',eye(neq)); kron(t.',eye(neq)) ];
v2 = [((-1).^(1:neq).');repmat(((-1).^(0:neq-1).'),d-1,1)];
M2 = diag(v2) * M1;
% M2 = diag([((-1).^(1:neq).');((-1).^(0:neq-1).')]) * M1;

% simplify( inv(v)*det(v) == [ cross( v(2,:), v(3,:) ).', ...
%                              cross( v(3,:), v(1,:) ).', ...
%                              cross( v(1,:), v(2,:) ).' ] )
% thus inv(v) == v.'

vec = @(x) x(:);

A1 = M2\M1;

e1 = zeros(d,1); e1(1) = 0;
% A2 = v.'*(eye(d) - 2*diag(e1))*v - diag(sum(v.*v,2));
% A2 = v.'*(eye(d) - 2*diag(e1))*v - eye(d);
A2 = @(p) v.'*diag( ones(d,1) - (1 +(-1)^p) * unit_vec(d,1) )*v - eye(d);

simplify( A2(2)-A2(1) == -2*n*(n.') )

% for i=1:3
%     I = circshift(eye(d),i-1,1);
%     A_1 = subs( A1, num2cell( vec(symvar(A1)) ), num2cell( I(:) ) );
%     A_2 = subs( A2, num2cell( vec(symvar(A2)) ), num2cell( I(:) ) );
%     all(A_1==A_2,"all")
% end


function e = unit_vec(d,n)
I = eye(d);
e = I(:,n);
end