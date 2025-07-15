%% (continued) prototype script for slip-wall BC for variational rec. (07/01/2025)
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
neq = 4;
syms n [d 1] real
syms u [d 1] real
syms P [d d] real
assume(dot(n,n)==1)
for j = 1:d
    for i = 1:d
        P(i,j) = n(i)*n(j);
    end
end
% P = P / (n.'*n);
% P = simplify(n*pinv(n));
R = eye(d) - P;

simplify(P*u + R*u)

% simplify( ((R-P)*u).*((R-P)*ones(3,1)) == ((R-P).*((R-P)*ones(3,1)))*u )

% simplify( (u - (R-P)*u).^2 == ( (eye(d) - (R-P))*u).^2 )

M1 = (eye(d) - (R-P));

M2 = (eye(d) + (R-P));


% simplify( (u - (R-P)*u).^2 == diag((M1*u)*((M1*u).') ) )

% simplify( (u + (R-P)*u).^2 == diag((M2*u)*((M2*u).') ) )

% M1_M1 = [ kron(M1(1,:),M1(1,:)); ...
%           kron(M1(2,:),M1(2,:)); ...
%           kron(M1(3,:),M1(3,:)) ];
% M2_M2 = [ kron(M2(1,:),M2(1,:)); ...
%           kron(M2(2,:),M2(2,:)); ...
%           kron(M2(3,:),M2(3,:)) ];



B1 = [ kron(P(1,:),P(1,:)); ...
       kron(P(2,:),P(2,:)); ...
       kron(P(3,:),P(3,:)) ];

B2 = split_face_square_n(n);

B3 = split_face_square_n_2(n);


M1_M1 = 4*[ kron(P(1,:),P(1,:)); ...
            kron(P(2,:),P(2,:)); ...
            kron(P(3,:),P(3,:)) ];
M2_M2 = 4*[ kron(R(1,:),R(1,:)); ...
            kron(R(2,:),R(2,:)); ...
            kron(R(3,:),R(3,:)) ];

simplify( (u - (R-P)*u).^2 == M1_M1*kron(u,u) )

simplify( (u + (R-P)*u).^2 == M2_M2*kron(u,u) )


% simplify( ([ kron(P(1,1:2),P(1,1:2)); kron(P(2,1:2),P(2,1:2)) ]*kron(u(1:2),u(1:2)) ) / ( n(1:2).' * u(1:2) ) )

% [ kron(P(1,1:2),P(1,1:2)); kron(P(2,1:2),P(2,1:2)) ]*[2*u1, 0; u2, u1; u2, u1; 0, 2*u2 ]

% 2*expand(diag(n(1:2))*( (n(1:2))*(n(1:2).') ) * (n(1:2).' * u(1:2) ))

% simplify( 2*expand(diag(n(1:2))*( (n(1:2))*(n(1:2).') ) * (n(1:2).' * u(1:2) )) == [ kron(P(1,1:2),P(1,1:2)); kron(P(2,1:2),P(2,1:2)) ]*[2*u1, 0; u2, u1; u2, u1; 0, 2*u2 ] )



function M = split_face_square_n(n)
d = numel(n);
syms M [d,d^2]

% for k = 1:d
%     for j = 1:d
%         for i = 1:d
%             M(k,(i-1)*d+j) = n(k)^2 * n(i) * n(j);
%         end
%     end
% end

for i = 1:d
    for j = 1:d
        for k = 1:d
            M(k,(i-1)*d+j) = n(k)^2 * n(i) * n(j);
        end
    end
end

end

function M = split_face_square_n_2(n)
M = kron(n.',n*n.').*n;
end

function D = delta_addition(d)
% delta_ik + delta_jk
D = zeros(d^2,d);
for k = 1:d
    for i = 1:d
        for j = 1:d
            if (i==k)
                D((i-1)*d+j,k) = D((i-1)*d+j,k) + 1;
            end
            % if (j==k)
            %     D((i-1)*d+j,k) = D((i-1)*d+j,k) + 1;
            % end
        end
    end
end
end