%% Using sparse matrix for generating LHS (06/23/2026)
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

% X = ones(10,10);
% 
% [r,c] = size(X);
% i     = 1:numel(X);
% j     = repmat(1:c,r,1);
% B     = sparse(i',j(:),X(:));
% 
% spy(B);

M = 10; n = 3;  p = 5;
% i     = 1:M*n;
A = zeros(M*n);
% (1:p)-(p+1)/2 + j

I = cell(p,1);
J = cell(p,1);
for k = 1:p
    d = k-(p+1)/2;
    [I{k},J{k}] = blk_diag_idx(M,n,d);
end


% [I1,J1] = blk_diag_idx(M,n,-2);
% [I2,J2] = blk_diag_idx(M,n,-1);
% [I3,J3] = blk_diag_idx(M,n, 0);
% I = [I1(:);I2(:);I3(:)];
% J = [J1(:);J2(:);J3(:)];
% A(I(:),I(:)) = 1;
A = sparse([I{:}],[J{:}],1);

A2 = create_sparse_banded_block_matrix(M,n,p,1);

% [I,J] = blk_diag_idx(M,n,+1);
% A(I(:),J(:)) = 1;


spy(A)



function [I,J] = blk_idx(n,i,d)
    i1 = n*(i-1)+1;
    i2 = i1 + n - 1;
    j1 = n*(i+d-1)+1;
    j2 = j1 + n - 1;
    [I,J] = ndgrid(i1:i2,j1:j2);
end

function [I,J] = blk_diag_idx(m,n,d)
    n_blk_rows     = m-abs(d);
    blk_row_offset = abs(min(d,0))+1;
    I = zeros(n,n,n_blk_rows);
    J = zeros(n,n,n_blk_rows);
    cnt = 0;
    for k = max(1,1-d):min(m,m-d)
        cnt = cnt + 1;
        i1 = n*(k-1)+1;
        i2 = i1 + n - 1;
        j1 = n*(k+d-1)+1;
        j2 = j1 + n - 1;
        % [I(:,:,k+1-blk_row_offset),J(:,:,k+1-blk_row_offset)] = ndgrid(i1:i2,j1:j2);
        [I(:,:,cnt),J(:,:,cnt)] = ndgrid(i1:i2,j1:j2);
    end
    I = I(:).';
    J = J(:).';
end

%% Set LHS Based on Grid Type
%% Apply Steger-Chaussee hyperbolic grid marching, eq 12
% Kinsey and Barth pg 5

% for i = 2:imax-1

% for i = 2:M-1
%     i1 = n*(i-1)+1:n*i;
%     j1 = n*(j-1) + n*(i-1)+1:n*i;
%   lhs( n*(i-1)+1:n*i, n*(i-1)+1:n*i 2*i-3:2*i-2) = 1;
%   lhs(2*i-1:2*i,2*i-3:2*i-2) = 1;
%   lhs(2*i-1:2*i,2*i+1:2*i+2) = 1;
% end
% % Lowest block diagonal
% for i = 3:imax-1
%   lhs(2*i-1,2*i-5) = lhs(2*i-1,2*i-5) + 1;
%   lhs(2*i,  2*i-4) = lhs(2*i,  2*i-4) + 1;
% end
% 
% % Lower block diagonal
% i = 2;
% lhs(2*i-1,2*i-3) = lhs(2*i-1,2*i-3) + 1;
% lhs(2*i,  2*i-2) = lhs(2*i,  2*i-2) + 1;
% for i = 3:imax-1
%   lhs(2*i-1,2*i-3) = lhs(2*i-1,2*i-3) + 1;
%   lhs(2*i,  2*i-2) = lhs(2*i,  2*i-2) + 1;
% end
% 
% % Diagonal
% i = 2;
% lhs(2*i-1,2*i-1) = lhs(2*i-1,2*i-1) + 1;
% lhs(2*i,  2*i)   = lhs(2*i,  2*i)   + 1;
% for i = 3:imax-2
%   lhs(2*i-1,2*i-1) = lhs(2*i-1,2*i-1) + 1;
%   lhs(2*i,  2*i)   = lhs(2*i,  2*i)   + 1;
% end
% i = imax-1;
% lhs(2*i-1,2*i-1) = lhs(2*i-1,2*i-1) + 1;
% lhs(2*i,  2*i)   = lhs(2*i,  2*i)   + 1;
% 
% % Upper block diagonal
% for i = 2:imax-2
%   lhs(2*i-1,2*i+1) = lhs(2*i-1,2*i+1) + 1;
%   lhs(2*i,  2*i+2) = lhs(2*i,  2*i+2) + 1;
% end
% i = imax-1;
% lhs(2*i-1,2*i+1) = lhs(2*i-1,2*i+1) + 1;
% lhs(2*i,  2*i+2) = lhs(2*i,  2*i+2) + 1;
% 
% % Last diagonal
% for i = 2:imax-2
%   lhs(2*i-1,2*i+3) = lhs(2*i-1,2*i+3) + 1;
%   lhs(2*i,  2*i+4) = lhs(2*i,  2*i+4) + 1;
% end
% 
% %% Boundary Condition
% lhs(2,4) = 1;
% lhs(2,6) = 1;
% lhs(2*imax,2*imax-2) = 1;
% lhs(2*imax,2*imax-4) = 1;
% 
% lhs = sparse(lhs);

% M = 100; n = 2;
% a = rand((M-2)*n, n); b = rand((M-1)*n, n); c = 5+rand(M*n, n);  
% d = rand((M-1)*n, n); e = rand((M-2)*n, n); f = rand(M*n,10);
% x1 = BlockPentSolve(a, b, c, d, e, f);

function x = BlockPentSolve(a, b, c, d, e, f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BlockPentSolve.m
%
% Solves the block pentadiagonal system Ax = f, where a, b, c, d, and e are 
% the five diagonals of A. if the size of the block is n-by-n, then c is M-by-n, 
% b and d are each (M-n)-by-n, while a and e are each (M - 2n)-by-n, f is 
% M-by-X, M is divisible by n, n >= 1. 
%
% M = 12; n = 3;
% a = rand((M-2)*n, n); b = rand((M-1)*n, n); c = 5+rand(M*n, n);  
% d = rand((M-1)*n, n); e = rand((M-2)*n, n); f = rand(M*n,10);
% x = BlockPentSolve(a, b, c, d, e, f);
%
% We can make block matrix A = BlockDiag(a, n, n, -2) + 
% BlockDiag(b, n, n, -1) + BlockDiag(c, n, n) + BlockDiag(d, n, n, 1) + 
% BlockDiag(e, n, n, 2); and then x = A\d to confirm the solution.
%
% Computational Cost of this method is 2M(5n^3+4n^2-n/3). The cost is
% of order M*n^3. This is better than the backslash which is order (M*n)^3
%
% For M = 500, n = 4; This function is faster than inbuilt backslash by
% factor of 2.5
%
% Written by: Lateef Adewale Kareem    05/25/2022
% Contact:    talk2laton@yahoo.co.uk
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(d, 2); M = size(f, 1)/n;
blk  = @(i) ((i-1)*n + 1):(i*n);
for i = 3:M
    w1 = a(blk(i-2), :)/c(blk(i-2), :);
    b(blk(i-1), :) = b(blk(i-1), :) - w1*d(blk(i-2), :);
    c(blk(i), :) = c(blk(i), :) - w1*e(blk(i-2), :);
    f(blk(i), :) = f(blk(i), :) - w1*f(blk(i-2), :);
    w2 = b(blk(i-2), :)/c(blk(i-2), :);
    c(blk(i-1), :) = c(blk(i-1), :) - w2*d(blk(i-2), :);
    d(blk(i-1), :) = d(blk(i-1), :) - w2*e(blk(i-2), :);
    f(blk(i-1), :) = f(blk(i-1), :) - w2*f(blk(i-2), :);
end
w2 = b(blk(M-1), :)/c(blk(M-1), :);
c(blk(M), :) = c(blk(M), :) - w2*d(blk(M-1), :);
f(blk(M), :) = f(blk(M), :) - w2*f(blk(M-1), :);
x = 0*f;  x(blk(M), :) = c(blk(M), :)\f(blk(M), :);
x(blk(M-1), :) = c(blk(M-1), :)\(f(blk(M-1), :) - ...
    d(blk(M-1), :)*x(blk(M), :));
for i = M-2:-1:1
    x(blk(i), :) = c(blk(i), :)\(f(blk(i), :) - ...
        d(blk(i), :)*x(blk(i+1), :) - e(blk(i), :)*x(blk(i+2), :));   
end

end