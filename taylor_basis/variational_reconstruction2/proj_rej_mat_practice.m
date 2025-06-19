%% projection/rejection matrix (06/12/2025)
clc; clear; close all;

A = [4;9;-2];
% A = eye(3);

vec = [5;4;0];

P = A*((A.'*A)\(A.'));


n = size(A,1);

mag_inv = 1/dot(A,A);

R = eye(n) - P;

% R2 = eye(n);
% for j = 1:n
%     for i = 1:n
%         R2(i,j) = R2(i,j) - A(i)*A(j)*mag_inv;
%     end
% end
R2 = zeros(n);
for i = 1:n
    R2(i,i) = 1;
    for j = 1:i
        R2(i,j) = R2(i,j) - A(i)*A(j)*mag_inv;
        R2(j,i) = R2(i,j);
    end
end

R*vec