%% Number of terms in polynomial refresher (09/28/2026)
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

d = 1;
k = 1;



function val = n_weak_compositions( n, k )
val = nchoosek( n + k - 1, n );
% alternatively
% val = nchoosek( n + k - 1, k - 1 );
end

function val = n_terms( d, k )
% total number of terms in a polynomial of d variables and total degree k
val = nchoosek( d + k, d );

% alternatively:
% val = nchoosek( d + k, k );
end

function val = n_terms_sum( d, k )
val = sum( arrayfun(@(k)n_weak_compositions(k,d),0:k ) );
end

% f1 = @(d,k) nchoosek(d+k,k)
% f3 = @(d,k) sum( arrayfun(@(k)f2(d,k),0:k ) )
% f2 = @(d,k) nchoosek(d+k-1,k)