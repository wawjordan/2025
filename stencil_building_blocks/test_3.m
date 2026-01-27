%% ========================================================================
%  JALALI'S CURVILINEAR K-EXACT RECONSTRUCTION
%  Clean Implementation: Polynomial Mapping + Cubic Spline Wall Distance
%  
%  Reference: Jalali (2017) PhD Dissertation, Section 3.3.2
% =========================================================================
clc; clear; close all;

%% ========== PATH SETUP ==========
parent_dir_str = 'stencil_building_blocks';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts

%% ========== CONTROL PARAMETERS ==========
VERBOSE_OUTPUT = true;
MAP_ORDER      = 3;       % Polynomial order for mapping (3 = cubic)
ERROR_NORM     = 'L2';    % Options: 'L1', 'L2', 'Linf'

%% ========== GRID LOADING ==========
load_file = true;
agglom    = false;
grid_file = 'svf.grd';
GRID = load_gen_grid_for_testing(parent_dir, grid_file, agglom, load_file);

%% ========== BUILD WALL SPLINES ==========
blk = 1;
wall_splines = build_wall_splines(GRID, blk);

%% ========== STENCIL PARAMETERS ==========
idx = [1, 10, 1];       % Reference cell [i, j, k]
n_stencil = 26;         % Number of stencil cells
balanced = true;
order = 4;              % Reconstruction order

block_info_list(1) = block_info_t();
block_info_list(1).block_id = 1;
block_info_list(1).Ncells(:) = GRID.gblock(1).Ncells;

%% ========== BUILD STENCIL ==========
% [stencil_idx, n_stencil_out, stencil_cells] = cell_t.build_stencil_alt(...
%     blk, idx, n_stencil, block_info_list, balanced, false);
[stencil_idx, n_stencil_out, stencil_cells] = cell_t.build_stencil(...
    blk, idx, n_stencil, GRID, block_info_list, balanced, false);

%% ========== COMPUTE WALL DISTANCE FOR REFERENCE CELL ==========
[D_i, t_hat, n_hat, closest_pt] = compute_wall_distance(GRID, blk, idx, wall_splines);

%% ========== COMPUTE CURVILINEAR MOMENTS WITH POLYNOMIAL MAPPING ==========
[xbar_tn, xhat_tn, t_centroids, n_centroids, mapping_data] = ...
    compute_2D_moments_curvilinear_poly(GRID, blk, idx, stencil_idx, order, wall_splines, MAP_ORDER);

%% ========== ASSEMBLE AND SOLVE LS SYSTEM ==========
[A, b, col_scale] = assemble_tn_ls_system(GRID, stencil_idx, xhat_tn, order);
coeffs_scaled = A \ b;
coeffs = coeffs_scaled ./ col_scale';

%% ========== PRINT RESULTS ==========
fprintf('\n========== CURVILINEAR K-EXACT RECONSTRUCTION ==========\n\n');
fprintf('Reference cell: blk=%d, idx=[%d,%d,%d]\n', blk, idx(1), idx(2), idx(3));
fprintf('Wall distance D_i = %.6e\n', D_i);
fprintf('Wall tangent t_hat = [%.6f, %.6f]\n', t_hat(1), t_hat(2));
fprintf('Wall normal  n_hat = [%.6f, %.6f]\n', n_hat(1), n_hat(2));
fprintf('Closest wall point = [%.6f, %.6f]\n\n', closest_pt(1), closest_pt(2));
fprintf('Method: POLYNOMIAL MAPPING + CUBIC SPLINE\n');
fprintf('Map order: %d\n\n', MAP_ORDER);

% fprintf('--- Self moments in (t,n) frame ---\n');
% fprintf('xbar_tn(1,1) = %.10f  (should be 1.0)\n', xbar_tn(1,1));
% fprintf('xbar_tn(2,1) = %.6e  (t-moment, should be ~0)\n', xbar_tn(2,1));
% fprintf('xbar_tn(1,2) = %.6e  (n-moment, should be ~0)\n\n', xbar_tn(1,2));

if VERBOSE_OUTPUT
    fprintf('--- Stencil (t,n) centroid coordinates ---\n');
    fprintf('  Cell   t_centroid      n_centroid\n');
    for s = 1:n_stencil_out
        fprintf('  %2d     %+.6e   %+.6e\n', s, t_centroids(s), n_centroids(s));
    end
    fprintf('\n');
end

fprintf('--- LS System ---\n');
fprintf('Size: %d x %d\n', size(A,1), size(A,2));
fprintf('Condition number (scaled): %.4e\n', cond(A));
fprintf('Rank: %d / %d\n\n', rank(A), min(size(A)));

if VERBOSE_OUTPUT
    fprintf('--- Reconstruction coefficients ---\n');
    coeff_idx = 1;
    for n = 0:order
        for m = 0:(order-n)
            fprintf('a(%d,%d) = %+.6e\n', n, m, coeffs(coeff_idx));
            coeff_idx = coeff_idx + 1;
        end
    end
    fprintf('\n');
end

%% ========== VERIFY K-EXACTNESS (j=1 row, each cell gets own reconstruction) ==========
fprintf('--- Verification: L2 error at j=1 row (each cell has own reconstruction) ---\n');

imax = GRID.gblock(blk).Ncells(1);
j_level = 1;

if VERBOSE_OUTPUT
    fprintf('  i      Exact (quad)    Cart L2 Err     Curv L2 Err\n');
end

max_err_cart = 0;
max_err_curv = 0;
total_error_sq_cart = 0;
total_error_sq_curv = 0;
total_volume = 0;

for i = 1:imax
    idx_i = [i, j_level, 1];
    
    % Build stencil for THIS cell
    % [stencil_idx_i, ~, ~] = cell_t.build_stencil_alt(blk, idx_i, n_stencil, block_info_list, balanced, false);
    [stencil_idx_i, ~, ~] = cell_t.build_stencil(blk, idx_i, n_stencil, GRID, block_info_list, balanced, false);
    
    % Cartesian reconstruction for THIS cell
    [~, xhat_cart_i] = compute_2D_moments_cartesian(GRID, blk, idx_i, stencil_idx_i, order);
    [A_cart_i, b_cart_i, col_scale_cart_i] = assemble_cartesian_ls_system(GRID, stencil_idx_i, xhat_cart_i, order);
    coeffs_cart_i = (A_cart_i \ b_cart_i) ./ col_scale_cart_i';
    
    % Curvilinear reconstruction for THIS cell
    [~, xhat_tn_i, ~, ~, mapping_data_i] = compute_2D_moments_curvilinear_poly(GRID, blk, idx_i, stencil_idx_i, order, wall_splines, MAP_ORDER);
    [A_curv_i, b_curv_i, col_scale_curv_i] = assemble_tn_ls_system(GRID, stencil_idx_i, xhat_tn_i, order);
    coeffs_curv_i = (A_curv_i \ b_curv_i) ./ col_scale_curv_i';
    
    % Reference point for Cartesian
    x_i = GRID.gblock(blk).grid_vars.cell_c(1, idx_i(1), idx_i(2), idx_i(3));
    y_i = GRID.gblock(blk).grid_vars.cell_c(2, idx_i(1), idx_i(2), idx_i(3));
    
    % Evaluate error ONLY at this cell
    Q = GRID.gblock(blk).grid_vars.quad(idx_i(1), idx_i(2), idx_i(3));
    
    cell_error_sq_cart = 0;
    cell_error_sq_curv = 0;
    cell_volume = 0;
    exact_integral = 0;
    
    for q = 1:Q.n_quad
        x_q = Q.quad_pts(1, q);
        y_q = Q.quad_pts(2, q);
        w_q = Q.quad_wts(q);
        
        % Exact value
        f_exact = test_function(x_q, y_q);
        
        % Cartesian reconstruction
        dx = x_q - x_i;
        dy = y_q - y_i;
        f_recon_cart = evaluate_polynomial(coeffs_cart_i, dx, dy, order);
        
        % Curvilinear reconstruction
        [t_q, n_q] = evaluate_curvilinear_mapping(x_q, y_q, mapping_data_i);
        f_recon_curv = evaluate_polynomial(coeffs_curv_i, t_q, n_q, order);
        
        % Accumulate
        cell_error_sq_cart = cell_error_sq_cart + (f_exact - f_recon_cart)^2 * w_q;
        cell_error_sq_curv = cell_error_sq_curv + (f_exact - f_recon_curv)^2 * w_q;
        cell_volume = cell_volume + w_q;
        exact_integral = exact_integral + f_exact * w_q;
    end
    
    exact_avg = exact_integral / cell_volume;
    cell_L2_err_cart = sqrt(cell_error_sq_cart / cell_volume);
    cell_L2_err_curv = sqrt(cell_error_sq_curv / cell_volume);
    
    max_err_cart = max(max_err_cart, cell_L2_err_cart);
    max_err_curv = max(max_err_curv, cell_L2_err_curv);
    
    total_error_sq_cart = total_error_sq_cart + cell_error_sq_cart;
    total_error_sq_curv = total_error_sq_curv + cell_error_sq_curv;
    total_volume = total_volume + cell_volume;
    
    if VERBOSE_OUTPUT
        fprintf('  %2d     %+.6e   %.4e        %.4e\n', i, exact_avg, cell_L2_err_cart, cell_L2_err_curv);
    end
end

global_L2_err_cart = sqrt(total_error_sq_cart / total_volume);
global_L2_err_curv = sqrt(total_error_sq_curv / total_volume);

fprintf('\n                        Cartesian       Curvilinear\n');
fprintf('Max cell L2 error:      %.4e      %.4e\n', max_err_cart, max_err_curv);
fprintf('Global L2 error:        %.4e      %.4e\n', global_L2_err_cart, global_L2_err_curv);

if max_err_cart < 1e-10
    fprintf('Cartesian:   SUCCESS (k-exactness verified)\n');
end
if max_err_curv < 1e-10
    fprintf('Curvilinear: SUCCESS (k-exactness verified)\n');
end

%% ========== COMPARE CONDITION NUMBERS ==========
fprintf('\n--- Condition Number Comparison ---\n');
[kappa_cart, kappa_cart_scaled, kappa_curv, kappa_curv_scaled] = ...
    compare_condition_numbers(GRID, blk, idx, stencil_idx, order, wall_splines, MAP_ORDER);

fprintf('Cartesian (unscaled):   %.4e\n', kappa_cart);
fprintf('Cartesian (scaled):     %.4e\n', kappa_cart_scaled);
fprintf('Curvilinear (unscaled): %.4e\n', kappa_curv);
fprintf('Curvilinear (scaled):   %.4e\n', kappa_curv_scaled);
fprintf('Improvement factor:     %.1fx\n', kappa_cart_scaled / kappa_curv_scaled);

%% ========== VISUALIZATION ==========
% Plot stencil
figure;
plot_stencil_cells_2D(GRID, stencil_idx, 'k');
label_stencil_cells_degree(GRID, stencil_cells(1:n_stencil_out));
title('Stencil cells');
grid on;

% Plot coordinate systems
visualize_coordinate_systems(GRID, blk, idx, stencil_idx, t_hat, n_hat, wall_splines);

%% ========== MULTI-GRID ANALYSIS ==========
grid_files = { 'svf33.grd', 'svf65.grd', 'svf129.grd', 'svf257.grd'};
%grid_files = {'flat_plate.35x25.grd','flat_plate.69x49.grd', 'flat_plate.137x97.grd','flat_plate.273x193.grd'};
plot_condition_number_vs_h(grid_files, blk, 1, n_stencil, order, MAP_ORDER);
% plot_error_norm_vs_h(grid_files, blk, 1, n_stencil, order, MAP_ORDER, ERROR_NORM);
% plot_convergence_j1_row(grid_files, blk, n_stencil, order, MAP_ORDER);
plot_convergence_j11_row(grid_files, blk, n_stencil, order, MAP_ORDER);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ==================== GRID LOADING ====================
function GRID = load_gen_grid_for_testing(parent_dir, grid_file, agglom, load_file)
    grid_t_file = [extractBefore(grid_file, '.'), '_GRID'];
    grid_file_path = fullfile(parent_dir, grid_file);
    if agglom
        grid_t_file = fullfile(parent_dir, [grid_t_file, '_agglom.mat']);
    else
        grid_t_file = fullfile(parent_dir, [grid_t_file, '.mat']);
    end
    if isfile(grid_t_file) && load_file
        load(grid_t_file, "GRID");
    else
        GRID = read_grd_file_to_struct(grid_file_path);
        if agglom
            GRID = grid_type(GRID, agglomerate=true, calc_quads=true, nquad=4, nskip=[2,2,1]);
        else
            GRID = grid_type(GRID, calc_quads=true, nquad=4);
        end
        save(grid_t_file, "GRID");
    end
end

function GRID = load_grid_for_plotting(parent_dir, grid_file)
    grid_t_file = [extractBefore(grid_file, '.'), '_GRID'];
    grid_path = fullfile(parent_dir, grid_file);
    grid_t_file = fullfile(parent_dir, [grid_t_file, '.mat']);
    if isfile(grid_t_file)
        load(grid_t_file, "GRID");
    else
        GRID = read_grd_file_to_struct(grid_path);
        GRID = grid_type(GRID, calc_quads=true, nquad=4);
        save(grid_t_file, "GRID");
    end
end

%% ==================== CUBIC SPLINE CONSTRUCTION (Jalali Eq. 3.17) ====================
function wall_splines = build_wall_splines(GRID, blk)
    jmax = GRID.gblock(blk).jmax;
    wall_splines.j1 = build_single_wall_spline(GRID, blk, 1);
    wall_splines.jmax = build_single_wall_spline(GRID, blk, jmax);
end

function spline = build_single_wall_spline(GRID, blk, wall_j)
    imax = GRID.gblock(blk).imax;
    x_wall = zeros(imax, 1);
    y_wall = zeros(imax, 1);
    for i = 1:imax
        x_wall(i) = GRID.gblock(blk).x(i, wall_j, 1);
        y_wall(i) = GRID.gblock(blk).y(i, wall_j, 1);
    end
    
    spline.n_seg = imax - 1;
    spline.x_nodes = x_wall;
    spline.y_nodes = y_wall;
    spline.alpha = zeros(spline.n_seg, 4);
    spline.beta = zeros(spline.n_seg, 4);
    
    for i = 1:spline.n_seg
        [spline.alpha(i,:), spline.beta(i,:)] = compute_hermite_coefficients(x_wall, y_wall, i, imax);
    end
end

function [alpha, beta] = compute_hermite_coefficients(x, y, seg_idx, imax)
    x0 = x(seg_idx);
    x1 = x(seg_idx + 1);
    y0 = y(seg_idx);
    y1 = y(seg_idx + 1);
    
    if seg_idx == 1
        tx0 = x(2) - x(1);
        ty0 = y(2) - y(1);
    else
        tx0 = 0.5 * (x(seg_idx + 1) - x(seg_idx - 1));
        ty0 = 0.5 * (y(seg_idx + 1) - y(seg_idx - 1));
    end
    
    if seg_idx == imax - 1
        tx1 = x(imax) - x(imax - 1);
        ty1 = y(imax) - y(imax - 1);
    else
        tx1 = 0.5 * (x(seg_idx + 2) - x(seg_idx));
        ty1 = 0.5 * (y(seg_idx + 2) - y(seg_idx));
    end
    
    alpha = [2*x0 - 2*x1 + tx0 + tx1, ...
             -3*x0 + 3*x1 - 2*tx0 - tx1, ...
             tx0, x0];
    beta = [2*y0 - 2*y1 + ty0 + ty1, ...
            -3*y0 + 3*y1 - 2*ty0 - ty1, ...
            ty0, y0];
end

%% ==================== WALL DISTANCE (Jalali Eq. 3.18) ====================
function [D, t_hat, n_hat, closest_pt] = compute_wall_distance(GRID, blk, idx, wall_splines)
    x0 = GRID.gblock(blk).grid_vars.cell_c(1, idx(1), idx(2), idx(3));
    y0 = GRID.gblock(blk).grid_vars.cell_c(2, idx(1), idx(2), idx(3));
    
    [D1, t1, n1, pt1] = compute_wall_distance_spline(wall_splines.j1, x0, y0);
    [D2, t2, n2, pt2] = compute_wall_distance_spline(wall_splines.jmax, x0, y0);
    
    if D1 <= D2
        D = D1; t_hat = t1; n_hat = n1; closest_pt = pt1;
    else
        D = D2; t_hat = t2; n_hat = n2; closest_pt = pt2;
    end
end

function D = compute_wall_distance_scalar(x, y, wall_splines)
    [D1, ~, ~, ~] = compute_wall_distance_spline(wall_splines.j1, x, y);
    [D2, ~, ~, ~] = compute_wall_distance_spline(wall_splines.jmax, x, y);
    D = min(D1, D2);
    % D = D1;
end

function [D, t_hat, n_hat, closest_pt] = compute_wall_distance_spline(spline, x0, y0)
    candidate_segs = find_candidate_segments(spline, x0, y0);
    if isempty(candidate_segs)
        candidate_segs = 1:spline.n_seg;
    end
    
    min_dist_sq = inf;
    best_t = 0;
    best_seg = 1;
    
    for seg = candidate_segs
        [tc, dist_sq] = newton_closest_point(spline.alpha(seg,:), spline.beta(seg,:), x0, y0);
        if dist_sq < min_dist_sq
            min_dist_sq = dist_sq;
            best_t = tc;
            best_seg = seg;
        end
    end
    
    % Check endpoints
    for i = 1:spline.n_seg + 1
        dist_sq = (spline.x_nodes(i) - x0)^2 + (spline.y_nodes(i) - y0)^2;
        if dist_sq < min_dist_sq
            min_dist_sq = dist_sq;
            if i <= spline.n_seg
                best_seg = i; best_t = 0;
            else
                best_seg = spline.n_seg; best_t = 1;
            end
        end
    end
    
    D = sqrt(min_dist_sq);
    alpha = spline.alpha(best_seg, :);
    beta = spline.beta(best_seg, :);
    
    closest_pt = [alpha(1)*best_t^3 + alpha(2)*best_t^2 + alpha(3)*best_t + alpha(4); ...
                  beta(1)*best_t^3 + beta(2)*best_t^2 + beta(3)*best_t + beta(4)];
    
    dx_dt = 3*alpha(1)*best_t^2 + 2*alpha(2)*best_t + alpha(3);
    dy_dt = 3*beta(1)*best_t^2 + 2*beta(2)*best_t + beta(3);
    mag = sqrt(dx_dt^2 + dy_dt^2);
    
    if mag < 1e-30
        t_hat = [1; 0];
    else
        t_hat = [dx_dt; dy_dt] / mag;
    end
    
    n_hat = [-t_hat(2); t_hat(1)];
    if dot(n_hat, [x0; y0] - closest_pt) < 0
        n_hat = -n_hat;
    end
end

function candidate_segs = find_candidate_segments(spline, x0, y0)
    candidate_segs = [];
    for seg = 1:spline.n_seg
        alpha = spline.alpha(seg, :);
        beta = spline.beta(seg, :);
        
        % Normals at endpoints
        dx0 = alpha(3); dy0 = beta(3);
        mag0 = sqrt(dx0^2 + dy0^2);
        if mag0 < 1e-30, n1 = [0; 1]; else, n1 = [-dy0; dx0] / mag0; end
        
        dx1 = 3*alpha(1) + 2*alpha(2) + alpha(3);
        dy1 = 3*beta(1) + 2*beta(2) + beta(3);
        mag1 = sqrt(dx1^2 + dy1^2);
        if mag1 < 1e-30, n2 = [0; 1]; else, n2 = [-dy1; dx1] / mag1; end
        
        P1 = [spline.x_nodes(seg); spline.y_nodes(seg)];
        P2 = [spline.x_nodes(seg + 1); spline.y_nodes(seg + 1)];
        
        r1 = [x0; y0] - P1;
        r2 = [x0; y0] - P2;
        
        cross1 = r1(1)*n1(2) - r1(2)*n1(1);
        cross2 = r2(1)*n2(2) - r2(2)*n2(1);
        
        if cross1 * cross2 < 0
            candidate_segs = [candidate_segs, seg]; 
        end
    end
end

function [tc, dist_sq] = newton_closest_point(alpha, beta, x0, y0)
    t = 0.5;
    for iter = 1:50
        x_t = alpha(1)*t^3 + alpha(2)*t^2 + alpha(3)*t + alpha(4);
        y_t = beta(1)*t^3 + beta(2)*t^2 + beta(3)*t + beta(4);
        
        dx_dt = 3*alpha(1)*t^2 + 2*alpha(2)*t + alpha(3);
        dy_dt = 3*beta(1)*t^2 + 2*beta(2)*t + beta(3);
        
        d2x_dt2 = 6*alpha(1)*t + 2*alpha(2);
        d2y_dt2 = 6*beta(1)*t + 2*beta(2);
        
        f = (x_t - x0)*dx_dt + (y_t - y0)*dy_dt;
        df = dx_dt^2 + dy_dt^2 + (x_t - x0)*d2x_dt2 + (y_t - y0)*d2y_dt2;
        
        if abs(df) < 1e-30, break; end
        
        dt = -f / df;
        t = max(0, min(1, t + dt));
        
        if abs(dt) < 1e-14, break; end
    end
    
    tc = t;
    x_tc = alpha(1)*tc^3 + alpha(2)*tc^2 + alpha(3)*tc + alpha(4);
    y_tc = beta(1)*tc^3 + beta(2)*tc^2 + beta(3)*tc + beta(4);
    dist_sq = (x_tc - x0)^2 + (y_tc - y0)^2;
    
    % Check endpoints
    for t_test = [0, 1]
        x_t = alpha(1)*t_test^3 + alpha(2)*t_test^2 + alpha(3)*t_test + alpha(4);
        y_t = beta(1)*t_test^3 + beta(2)*t_test^2 + beta(3)*t_test + beta(4);
        d2 = (x_t - x0)^2 + (y_t - y0)^2;
        if d2 < dist_sq
            dist_sq = d2;
            tc = t_test;
        end
    end
end

%% ==================== PRINCIPAL AXES ====================
function [theta, R] = compute_principal_axes(xbar)
    M_xx = xbar(3, 1);
    M_xy = xbar(2, 2);
    M_yy = xbar(1, 3);
    M = [M_xx, M_xy; M_xy, M_yy];
    [V, D] = eig(M);
    eig_vals = diag(D);
    [~, sort_idx] = sort(eig_vals, 'descend');
    v1 = V(:, sort_idx(1));
    theta = atan2(v1(2), v1(1));
    c = cos(theta); s = sin(theta);
    R = [c, s; -s, c];
end

%% ==================== CARTESIAN MOMENTS ====================
function xbar = compute_cell_moments(GRID, blk, idx, order)
    xbar = zeros(order+1, order+1);
    Q = GRID.gblock(blk).grid_vars.quad(idx(1), idx(2), idx(3));
    xc = GRID.gblock(blk).grid_vars.cell_c(1, idx(1), idx(2), idx(3));
    yc = GRID.gblock(blk).grid_vars.cell_c(2, idx(1), idx(2), idx(3));
    
    quad_pts = Q.quad_pts;
    quad_wts = Q.quad_wts;
    n_quad = Q.n_quad;
    cell_area = sum(quad_wts);
    
    for n = 0:order
        for m = 0:(order - n)
            moment_sum = 0;
            for q = 1:n_quad
                dx = quad_pts(1, q) - xc;
                dy = quad_pts(2, q) - yc;
                moment_sum = moment_sum + (dx^n) * (dy^m) * quad_wts(q);
            end
            xbar(n+1, m+1) = moment_sum / cell_area;
        end
    end
end

function [xbar, xhat] = compute_2D_moments_cartesian(GRID, blk, idx, stencil_idx, order)
    n_stencil = size(stencil_idx, 2);
    xi = GRID.gblock(blk).grid_vars.cell_c(1, idx(1), idx(2), idx(3));
    yi = GRID.gblock(blk).grid_vars.cell_c(2, idx(1), idx(2), idx(3));
    
    xbar = compute_cell_moments(GRID, blk, idx, order);
    xhat = zeros(n_stencil, order+1, order+1);
    
    for s = 1:n_stencil
        blk_j = stencil_idx(1, s);
        idx_j = stencil_idx(2:4, s)';
        xj = GRID.gblock(blk_j).grid_vars.cell_c(1, idx_j(1), idx_j(2), idx_j(3));
        yj = GRID.gblock(blk_j).grid_vars.cell_c(2, idx_j(1), idx_j(2), idx_j(3));
        dx = xj - xi;
        dy = yj - yi;
        xbar_j = compute_cell_moments(GRID, blk_j, idx_j, order);
        
        for n = 0:order
            for m = 0:(order - n)
                tmp = 0;
                for k = 0:n
                    for l = 0:m
                        tmp = tmp + nchoosek(n,k) * nchoosek(m,l) * (dx^k) * (dy^l) * xbar_j(n-k+1, m-l+1);
                    end
                end
                xhat(s, n+1, m+1) = tmp;
            end
        end
    end
end

%% ==================== POLYNOMIAL MAPPING (Jalali Eq. 3.22-3.23) ====================
function [xbar_tn, xhat_tn, t_centroids, n_centroids, mapping_data] = ...
    compute_2D_moments_curvilinear_poly(GRID, blk_i, idx_i, stencil_idx, order, wall_splines, map_order)
    
    n_stencil = size(stencil_idx, 2);
    
    % Reference cell data
    x_i = GRID.gblock(blk_i).grid_vars.cell_c(1, idx_i(1), idx_i(2), idx_i(3));
    y_i = GRID.gblock(blk_i).grid_vars.cell_c(2, idx_i(1), idx_i(2), idx_i(3));
    D_i = compute_wall_distance_scalar(x_i, y_i, wall_splines);
    [~, t_hat_i, ~, ~] = compute_wall_distance(GRID, blk_i, idx_i, wall_splines);
    
    % Principal axes from Cartesian moments
    xbar_cart = compute_cell_moments(GRID, blk_i, idx_i, 2);
    [~, R] = compute_principal_axes(xbar_cart);
    
    % Compute polynomial mapping coefficients
    [t_coeffs, n_coeffs] = compute_mapping_coefficients(GRID, stencil_idx, ...
        wall_splines, map_order, R, t_hat_i, D_i, x_i, y_i);
    
    % Store mapping data for later use
    mapping_data.x_i = x_i;
    mapping_data.y_i = y_i;
    mapping_data.R = R;
    mapping_data.t_coeffs = t_coeffs;
    mapping_data.n_coeffs = n_coeffs;
    mapping_data.map_order = map_order;
    
    % Initialize outputs
    xhat_tn = zeros(n_stencil, order+1, order+1);
    t_centroids = zeros(n_stencil, 1);
    n_centroids = zeros(n_stencil, 1);
    
    % Self-moments for reference cell
    xbar_tn = compute_cell_self_moments_tn_poly(GRID, blk_i, idx_i, mapping_data, order);
    
    % Cross-cell moments via direct quadrature with polynomial mapping
    for s = 1:n_stencil
        blk_j = stencil_idx(1, s);
        idx_j = stencil_idx(2:4, s)';
        
        Q_j = GRID.gblock(blk_j).grid_vars.quad(idx_j(1), idx_j(2), idx_j(3));
        quad_pts_j = Q_j.quad_pts;
        quad_wts_j = Q_j.quad_wts;
        n_quad_j = Q_j.n_quad;
        cell_area_j = sum(quad_wts_j);
        
        % Centroid in (t,n)
        x_j = GRID.gblock(blk_j).grid_vars.cell_c(1, idx_j(1), idx_j(2), idx_j(3));
        y_j = GRID.gblock(blk_j).grid_vars.cell_c(2, idx_j(1), idx_j(2), idx_j(3));
        [t_cen_j, n_cen_j] = evaluate_curvilinear_mapping(x_j, y_j, mapping_data);
        
        t_centroids(s) = t_cen_j;
        n_centroids(s) = n_cen_j;
        
        % Direct quadrature
        moment_accum = zeros(order+1, order+1);
        for q = 1:n_quad_j
            x_q = quad_pts_j(1, q);
            y_q = quad_pts_j(2, q);
            w_q = quad_wts_j(q);
            
            [t_q, n_q] = evaluate_curvilinear_mapping(x_q, y_q, mapping_data);
            
            for n_pow = 0:order
                t_power = t_q^n_pow;
                for m_pow = 0:(order - n_pow)
                    moment_accum(n_pow+1, m_pow+1) = moment_accum(n_pow+1, m_pow+1) + ...
                        w_q * t_power * (n_q^m_pow);
                end
            end
        end
        
        xhat_tn(s, :, :) = moment_accum / cell_area_j;
    end
end

function [t_coeffs, n_coeffs] = compute_mapping_coefficients(GRID, stencil_idx, ...
    wall_splines, map_order, R, t_hat_i, D_i, x_i, y_i)
    
    n_stencil = size(stencil_idx, 2);
    n_terms = (map_order + 1) * (map_order + 2) / 2 - 1;  % Exclude constant
    
    A_aux = zeros(n_stencil, n_terms);
    t_rhs = zeros(n_stencil, 1);
    n_rhs = zeros(n_stencil, 1);
    
    for s = 1:n_stencil
        blk_j = stencil_idx(1, s);
        idx_j = stencil_idx(2:4, s)';
        
        x_j = GRID.gblock(blk_j).grid_vars.cell_c(1, idx_j(1), idx_j(2), idx_j(3));
        y_j = GRID.gblock(blk_j).grid_vars.cell_c(2, idx_j(1), idx_j(2), idx_j(3));
        D_j = compute_wall_distance_scalar(x_j, y_j, wall_splines);
        
        dx = x_j - x_i;
        dy = y_j - y_i;
        
        % Rotate to principal coordinates
        x_prime = R(1,1)*dx + R(1,2)*dy;
        y_prime = R(2,1)*dx + R(2,2)*dy;
        
        % Known curvilinear coordinates
        t_j = dx * t_hat_i(1) + dy * t_hat_i(2);
        n_j = D_j - D_i;
        
        % Fill matrix row 
        col = 1;
        for p = 1:map_order
            for q = 0:p
                A_aux(s, col) = (x_prime^(p-q)) * (y_prime^q);
                col = col + 1;
            end
        end
        
        t_rhs(s) = t_j;
        n_rhs(s) = n_j;
    end
    
    % Column scaling
    col_scale = max(abs(A_aux), [], 1);
    col_scale(col_scale < eps) = 1;
    A_aux_scaled = A_aux ./ col_scale;
    
    % Solve
    t_coeffs_scaled = A_aux_scaled \ t_rhs;
    n_coeffs_scaled = A_aux_scaled \ n_rhs;
    
    t_coeffs = t_coeffs_scaled ./ col_scale';
    n_coeffs = n_coeffs_scaled ./ col_scale';
end

function [t, n] = evaluate_curvilinear_mapping(x, y, mapping_data)
    dx = x - mapping_data.x_i;
    dy = y - mapping_data.y_i;
    R = mapping_data.R;
    
    x_prime = R(1,1)*dx + R(1,2)*dy;
    y_prime = R(2,1)*dx + R(2,2)*dy;
    
    t = 0; n = 0;
    col = 1;
    for p = 1:mapping_data.map_order
        for q = 0:p
            term = (x_prime^(p-q)) * (y_prime^q);
            t = t + mapping_data.t_coeffs(col) * term;
            n = n + mapping_data.n_coeffs(col) * term;
            col = col + 1;
        end
    end
end

function xbar_tn = compute_cell_self_moments_tn_poly(GRID, blk, idx, mapping_data, order)
    xbar_tn = zeros(order+1, order+1);
    
    Q = GRID.gblock(blk).grid_vars.quad(idx(1), idx(2), idx(3));
    x_self = GRID.gblock(blk).grid_vars.cell_c(1, idx(1), idx(2), idx(3));
    y_self = GRID.gblock(blk).grid_vars.cell_c(2, idx(1), idx(2), idx(3));
    
    quad_pts = Q.quad_pts;
    quad_wts = Q.quad_wts;
    n_quad = Q.n_quad;
    cell_area = sum(quad_wts);
    
    [t_self, n_self] = evaluate_curvilinear_mapping(x_self, y_self, mapping_data);
    
    for q = 1:n_quad
        x_q = quad_pts(1, q);
        y_q = quad_pts(2, q);
        w_q = quad_wts(q);
        
        [t_q, n_q] = evaluate_curvilinear_mapping(x_q, y_q, mapping_data);
        
        dt = t_q - t_self;
        dn = n_q - n_self;
        
        for n_pow = 0:order
            for m_pow = 0:(order - n_pow)
                xbar_tn(n_pow+1, m_pow+1) = xbar_tn(n_pow+1, m_pow+1) + ...
                    w_q * (dt^n_pow) * (dn^m_pow);
            end
        end
    end
    
    xbar_tn = xbar_tn / cell_area;
end

%% ==================== LS SYSTEM ASSEMBLY ====================
function [A, b, col_scale] = assemble_tn_ls_system(GRID, stencil_idx, xhat_tn, order)
    n_stencil = size(stencil_idx, 2);
    n_coeffs = (order + 1) * (order + 2) / 2;
    
    A = zeros(n_stencil, n_coeffs);
    b = zeros(n_stencil, 1);
    
    for s = 1:n_stencil
        blk_j = stencil_idx(1, s);
        idx_j = stencil_idx(2:4, s)';
        
        col = 1;
        for n = 0:order
            for m = 0:(order - n)
                A(s, col) = xhat_tn(s, n+1, m+1);
                col = col + 1;
            end
        end
        
        b(s) = compute_cell_average(GRID, blk_j, idx_j);
    end
    
    col_scale = max(abs(A), [], 1);
    col_scale(col_scale < eps) = 1;
    A = A ./ col_scale;
end

function [A, b, col_scale] = assemble_cartesian_ls_system(GRID, stencil_idx, xhat_cart, order)
    n_stencil = size(stencil_idx, 2);
    n_coeffs = (order + 1) * (order + 2) / 2;
    
    A = zeros(n_stencil, n_coeffs);
    b = zeros(n_stencil, 1);
    
    for s = 1:n_stencil
        blk_j = stencil_idx(1, s);
        idx_j = stencil_idx(2:4, s)';
        
        col = 1;
        for n = 0:order
            for m = 0:(order - n)
                A(s, col) = xhat_cart(s, n+1, m+1);
                col = col + 1;
            end
        end
        
        b(s) = compute_cell_average(GRID, blk_j, idx_j);
    end
    
    col_scale = max(abs(A), [], 1);
    col_scale(col_scale < eps) = 1;
    A = A ./ col_scale;
end

%% ==================== TEST FUNCTION AND CELL AVERAGE ====================
function avg = compute_cell_average(GRID, blk, idx)
    Q = GRID.gblock(blk).grid_vars.quad(idx(1), idx(2), idx(3));
    integral = 0;
    total_weight = 0;
    for q = 1:Q.n_quad
        x_q = Q.quad_pts(1, q);
        y_q = Q.quad_pts(2, q);
        w = Q.quad_wts(q);
        f_q = test_function(x_q, y_q);
        integral = integral + f_q * w;
        total_weight = total_weight + w;
    end
    avg = integral / total_weight;
end

function f = test_function(x, y)
    f = sin(2*pi*x)*cos(2*pi*y);
end

%% ==================== RECONSTRUCTION EVALUATION ====================
function avg = evaluate_reconstruction_average_poly(GRID, blk_j, idx_j, coeffs, mapping_data, order)
    Q = GRID.gblock(blk_j).grid_vars.quad(idx_j(1), idx_j(2), idx_j(3));
    
    integral = 0;
    total_weight = 0;
    
    for q = 1:Q.n_quad
        x_q = Q.quad_pts(1, q);
        y_q = Q.quad_pts(2, q);
        w_q = Q.quad_wts(q);
        
        [t_q, n_q] = evaluate_curvilinear_mapping(x_q, y_q, mapping_data);
        
        phi_q = evaluate_polynomial(coeffs, t_q, n_q, order);
        integral = integral + phi_q * w_q;
        total_weight = total_weight + w_q;
    end
    avg = integral / total_weight;
end

function avg = evaluate_reconstruction_average_cartesian(GRID, blk_j, idx_j, blk_i, idx_i, coeffs, order)
    Q = GRID.gblock(blk_j).grid_vars.quad(idx_j(1), idx_j(2), idx_j(3));
    
    x_i = GRID.gblock(blk_i).grid_vars.cell_c(1, idx_i(1), idx_i(2), idx_i(3));
    y_i = GRID.gblock(blk_i).grid_vars.cell_c(2, idx_i(1), idx_i(2), idx_i(3));
    
    integral = 0;
    total_weight = 0;
    
    for q = 1:Q.n_quad
        x_q = Q.quad_pts(1, q);
        y_q = Q.quad_pts(2, q);
        w_q = Q.quad_wts(q);
        
        dx = x_q - x_i;
        dy = y_q - y_i;
        
        phi_q = evaluate_polynomial(coeffs, dx, dy, order);
        integral = integral + phi_q * w_q;
        total_weight = total_weight + w_q;
    end
    avg = integral / total_weight;
end

function phi = evaluate_polynomial(coeffs, t, n, order)
    phi = 0;
    coeff_idx = 1;
    for p = 0:order
        for m = 0:(order - p)
            phi = phi + coeffs(coeff_idx) * (t^p) * (n^m);
            coeff_idx = coeff_idx + 1;
        end
    end
end

%% ==================== CONDITION NUMBER COMPARISON ====================
function [kappa_cart, kappa_cart_scaled, kappa_curv, kappa_curv_scaled] = ...
    compare_condition_numbers(GRID, blk, idx, stencil_idx, order, wall_splines, map_order)
    
    % Cartesian
    [~, xhat_cart] = compute_2D_moments_cartesian(GRID, blk, idx, stencil_idx, order);
    [A_cart, ~, col_scale_cart] = assemble_cartesian_ls_system(GRID, stencil_idx, xhat_cart, order);
    
    kappa_cart = cond(A_cart .* col_scale_cart);
    kappa_cart_scaled = cond(A_cart);
    
    % Curvilinear
    [~, xhat_tn, ~, ~, ~] = compute_2D_moments_curvilinear_poly(GRID, blk, idx, stencil_idx, order, wall_splines, map_order);
    [A_curv, ~, col_scale_curv] = assemble_tn_ls_system(GRID, stencil_idx, xhat_tn, order);
    
    kappa_curv = cond(A_curv .* col_scale_curv);
    kappa_curv_scaled = cond(A_curv);
end

%% ==================== PLOTTING FUNCTIONS ====================
function grid_in_stencil = get_stencil_grid(Npts, stencil, grid)
    n_stencil = size(stencil, 2);
    grid_in_stencil = cell(n_stencil, 3);
    for si = 1:n_stencil
        b = stencil(1, si);
        cell_idx = stencil(2:4, si)';
        X = get_local_interp_grid(grid.gblock(b).grid_vars.quad(cell_idx(1), cell_idx(2), cell_idx(3)), grid.gblock(b).dim, Npts);
        grid_in_stencil{si, 1} = X{1};
        grid_in_stencil{si, 2} = X{2};
        grid_in_stencil{si, 3} = X{3};
    end
end

function X = get_local_interp_grid(Q, n_dim, npts)
    n_quad = repmat(round(Q.n_quad.^(1/n_dim)), 1, n_dim);
    quad_ref = quad_t.create_quad_ref_2D(n_quad(1));
    n_plot_pts = ones(3, 1);
    n_plot_pts(1:n_dim) = npts;
    xtmp = reshape(Q.quad_pts(1, :), n_quad);
    ytmp = reshape(Q.quad_pts(2, :), n_quad);
    ztmp = reshape(Q.quad_pts(2, :), n_quad);
    L = lagrange_interpolant(n_quad(1) + 1);
    [x1, y1, z1] = L.map_grid_from_quad(xtmp, ytmp, ztmp, n_plot_pts(1), n_plot_pts(2), n_plot_pts(3), quad_ref);
    X = {x1, y1, z1};
end

function plot_stencil_cells_2D(GRID, stencil, varargin)
    Npts = 2;
    if any(GRID.nskip > 1)
        Npts = 21;
    end
    grid_in_stencil = get_stencil_grid(Npts, stencil, GRID);
    n_stencil = size(stencil, 2);
    hold on;
    for n = 1:n_stencil
        X = grid_in_stencil{n, 1}(:, :, 1);
        Y = grid_in_stencil{n, 2}(:, :, 1);
        Z = grid_in_stencil{n, 3}(:, :, 1);
        plot_2D_mesh_edge(X, Y, Z, varargin{:});
    end
end

function label_stencil_cells_degree(GRID, stencil_cells)
    n_stencil = numel(stencil_cells);
    for n = 1:n_stencil
        b = stencil_cells(n).idx(1);
        cell_idx = stencil_cells(n).idx(2:4);
        degree = stencil_cells(n).degree;
        label_cell(GRID.gblock(b).grid_vars.cell_c(:, cell_idx(1), cell_idx(2), cell_idx(3)), '%d(%d)', [n; degree], 'Color', 'r');
    end
end

function plot_2D_mesh_edge(X, Y, Z, varargin)
    hold on;
    Xp = [X(1:end,1); X(end,1:end)'; X(end:-1:1,end); X(1,end:-1:1)'];
    Yp = [Y(1:end,1); Y(end,1:end)'; Y(end:-1:1,end); Y(1,end:-1:1)'];
    Zp = [Z(1:end,1); Z(end,1:end)'; Z(end:-1:1,end); Z(1,end:-1:1)'];
    plot3(Xp, Yp, Zp, varargin{:});
end

function label_cell(X, fmt, val, varargin)
    text(X(1), X(2), X(3), sprintf(fmt, val), varargin{:});
end

function visualize_coordinate_systems(GRID, blk, idx, stencil_idx, t_hat, n_hat, wall_splines)
    n_stencil = size(stencil_idx, 2);
    
    x_i = GRID.gblock(blk).grid_vars.cell_c(1, idx(1), idx(2), idx(3));
    y_i = GRID.gblock(blk).grid_vars.cell_c(2, idx(1), idx(2), idx(3));
    D_i = compute_wall_distance_scalar(x_i, y_i, wall_splines);
    
    x_coords = zeros(n_stencil, 1);
    y_coords = zeros(n_stencil, 1);
    t_coords = zeros(n_stencil, 1);
    n_coords = zeros(n_stencil, 1);
    
    for s = 1:n_stencil
        blk_j = stencil_idx(1, s);
        idx_j = stencil_idx(2:4, s)';
        x_j = GRID.gblock(blk_j).grid_vars.cell_c(1, idx_j(1), idx_j(2), idx_j(3));
        y_j = GRID.gblock(blk_j).grid_vars.cell_c(2, idx_j(1), idx_j(2), idx_j(3));
        D_j = compute_wall_distance_scalar(x_j, y_j, wall_splines);
        
        x_coords(s) = x_j - x_i;
        y_coords(s) = y_j - y_i;
        t_coords(s) = (x_j - x_i) * t_hat(1) + (y_j - y_i) * t_hat(2);
        n_coords(s) = D_j - D_i;
    end
    
    figure('Position', [100, 100, 1200, 500]);
    
    subplot(1, 2, 1);
    hold on;
    plot(x_coords, y_coords, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    plot(0, 0, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
    for s = 1:n_stencil
        text(x_coords(s), y_coords(s), sprintf(' %d', s), 'FontSize', 10);
    end
    grid on; axis equal;
    xlabel('x - x_i'); ylabel('y - y_i');
    title('Cartesian coordinates');
    
    subplot(1, 2, 2);
    hold on;
    plot(t_coords, n_coords, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot(0, 0, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
    for s = 1:n_stencil
        text(t_coords(s), n_coords(s), sprintf(' %d', s), 'FontSize', 10);
    end
    grid on; axis equal;
    xlabel('t (tangential)'); ylabel('n (wall-normal)');
    title('Curvilinear (t,n) - Polynomial Mapping');
    
    sgtitle(sprintf('Stencil coordinates: cell [%d,%d,%d]', idx(1), idx(2), idx(3)));
end

%% ====================  CONDITION NUMBER PLOT ====================
function plot_condition_number_vs_h(grid_files, blk, j_level, n_stencil, order, map_order)
    n_grids = length(grid_files);
    
    h_values = zeros(n_grids, 1);
    kappa_cart = zeros(n_grids, 1);
    kappa_cart_scaled = zeros(n_grids, 1);
    kappa_curv = zeros(n_grids, 1);
    kappa_curv_scaled = zeros(n_grids, 1);
    
    parent_dir = pwd;
    
    for g = 1:n_grids
        fprintf('Condition number - Processing grid %d/%d: %s\n', g, n_grids, grid_files{g});
        
        GRID = load_grid_for_plotting(parent_dir, grid_files{g});
        wall_splines = build_wall_splines(GRID, blk);
        
        i_mid = round(GRID.gblock(blk).Ncells(1) / 2);
        idx = [i_mid, j_level, 1];
        
        h_values(g) = compute_mesh_size(GRID, blk, idx);
        
        block_info_list(1) = block_info_t();
        block_info_list(1).block_id = 1;
        block_info_list(1).Ncells(:) = GRID.gblock(1).Ncells;
        
        [stencil_idx, ~, ~] = cell_t.build_stencil_alt(blk, idx, n_stencil, block_info_list, true, false);
        
        [kappa_cart(g), kappa_cart_scaled(g), kappa_curv(g), kappa_curv_scaled(g)] = ...
            compare_condition_numbers(GRID, blk, idx, stencil_idx, order, wall_splines, map_order);
    end
    
    % Sort by h (descending)
    [h_values, sort_idx] = sort(h_values, 'descend');
    kappa_cart = kappa_cart(sort_idx);
    kappa_cart_scaled = kappa_cart_scaled(sort_idx);
    kappa_curv = kappa_curv(sort_idx);
    kappa_curv_scaled = kappa_curv_scaled(sort_idx);
    
    % Plot
    figure('Position', [100, 100, 500, 400]);
    p1 = loglog(h_values, kappa_cart, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    hold on;
    p2 = loglog(h_values, kappa_cart_scaled, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    p3 = loglog(h_values, kappa_curv, 'b--s', 'LineWidth', 1.5, 'MarkerSize', 8);
    p4 = loglog(h_values, kappa_curv_scaled, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    
    % Format axes
    ax = gca;
    all_kappa = [kappa_cart; kappa_cart_scaled; kappa_curv; kappa_curv_scaled];
    y_min = floor(log10(min(all_kappa)));
    y_max = ceil(log10(max(all_kappa)));
    y_powers = y_min:2:y_max;
    ylim([10^y_min, 10^y_max]);
    yticks(10.^y_powers);
    yticklabels(arrayfun(@(p) sprintf('10^{%d}', p), y_powers, 'UniformOutput', false));
    set(ax, 'TickLabelInterpreter', 'tex');
    
    grid on; box on;
    xlabel('h', 'FontSize', 12);
    ylabel('Condition Number', 'FontSize', 12);
    title(sprintf('Condition Number vs Mesh Size (Order %d)', order), 'FontSize', 12);
    legend([p1, p2, p3, p4], {'Cartesian (unscaled)', 'Cartesian (scaled)', ...
        'Curvilinear (unscaled)', 'Curvilinear (scaled)'}, 'Location', 'best', 'FontSize', 9);
    
    % Print results
    fprintf('\n=== CONDITION NUMBER RESULTS ===\n');
    fprintf('%-15s %12s %12s %12s %12s %12s\n', 'Grid', 'h', 'Cart', 'Cart(sc)', 'Curv', 'Curv(sc)');
    fprintf('%s\n', repmat('-', 75, 1));
    for g = 1:n_grids
        fprintf('%-15s %12.4e %12.2e %12.2e %12.2e %12.2e\n', ...
            grid_files{sort_idx(g)}, h_values(g), kappa_cart(g), kappa_cart_scaled(g), kappa_curv(g), kappa_curv_scaled(g));
    end
end

%% ==================== ERROR NORM PLOT ====================
function plot_error_norm_vs_h(grid_files, blk, j_level, n_stencil, order, map_order, error_norm)
    n_grids = length(grid_files);
    
    h_values = zeros(n_grids, 1);
    err_cart = zeros(n_grids, 1);
    err_curv = zeros(n_grids, 1);
    err_curv_scaled = zeros(n_grids, 1);
    
    parent_dir = pwd;
    
    for g = 1:n_grids
        fprintf('Error norm - Processing grid %d/%d: %s\n', g, n_grids, grid_files{g});
        
        GRID = load_grid_for_plotting(parent_dir, grid_files{g});
        wall_splines = build_wall_splines(GRID, blk);
        
        i_mid = round(GRID.gblock(blk).Ncells(1) / 2);
        idx = [i_mid, j_level, 1];
        
        h_values(g) = compute_mesh_size(GRID, blk, idx);
        
        block_info_list(1) = block_info_t();
        block_info_list(1).block_id = 1;
        block_info_list(1).Ncells(:) = GRID.gblock(1).Ncells;
        
        [stencil_idx, n_stencil_out, ~] = cell_t.build_stencil_alt(blk, idx, n_stencil, block_info_list, true, false);
        
        % Cartesian reconstruction
        [~, xhat_cart] = compute_2D_moments_cartesian(GRID, blk, idx, stencil_idx, order);
        [A_cart, b_cart, col_scale_cart] = assemble_cartesian_ls_system(GRID, stencil_idx, xhat_cart, order);
        coeffs_cart_scaled = A_cart \ b_cart;
        coeffs_cart = coeffs_cart_scaled ./ col_scale_cart';
        
        % Curvilinear reconstruction
        [~, xhat_tn, ~, ~, mapping_data] = compute_2D_moments_curvilinear_poly(GRID, blk, idx, stencil_idx, order, wall_splines, map_order);
        [A_curv, b_curv, col_scale_curv] = assemble_tn_ls_system(GRID, stencil_idx, xhat_tn, order);
        
        % Unscaled solve
        A_curv_unscaled = A_curv .* col_scale_curv;
        coeffs_curv = A_curv_unscaled \ b_curv;
        
        % Scaled solve
        coeffs_curv_scaled_tmp = A_curv \ b_curv;
        coeffs_curv_scaled = coeffs_curv_scaled_tmp ./ col_scale_curv';
        
        % Compute errors
        errors_cart = zeros(n_stencil_out, 1);
        errors_curv = zeros(n_stencil_out, 1);
        errors_curv_scaled = zeros(n_stencil_out, 1);
        
        for s = 1:n_stencil_out
            blk_j = stencil_idx(1, s);
            idx_j = stencil_idx(2:4, s)';
            exact_avg = compute_cell_average(GRID, blk_j, idx_j);
            
            recon_cart = evaluate_reconstruction_average_cartesian(GRID, blk_j, idx_j, blk, idx, coeffs_cart, order);
            errors_cart(s) = abs(exact_avg - recon_cart);
            
            recon_curv = evaluate_reconstruction_average_poly(GRID, blk_j, idx_j, coeffs_curv, mapping_data, order);
            errors_curv(s) = abs(exact_avg - recon_curv);
            
            recon_curv_scaled = evaluate_reconstruction_average_poly(GRID, blk_j, idx_j, coeffs_curv_scaled, mapping_data, order);
            errors_curv_scaled(s) = abs(exact_avg - recon_curv_scaled);
        end
        
        % Compute norm
        switch error_norm
            case 'L1'
                err_cart(g) = mean(errors_cart);
                err_curv(g) = mean(errors_curv);
                err_curv_scaled(g) = mean(errors_curv_scaled);
            case 'L2'
                err_cart(g) = sqrt(mean(errors_cart.^2));
                err_curv(g) = sqrt(mean(errors_curv.^2));
                err_curv_scaled(g) = sqrt(mean(errors_curv_scaled.^2));
            case 'Linf'
                err_cart(g) = max(errors_cart);
                err_curv(g) = max(errors_curv);
                err_curv_scaled(g) = max(errors_curv_scaled);
        end
    end
    
    % Sort by h (descending)
    [h_values, sort_idx] = sort(h_values, 'descend');
    err_cart = err_cart(sort_idx);
    err_curv = err_curv(sort_idx);
    err_curv_scaled = err_curv_scaled(sort_idx);
    
    % Plot
    figure('Position', [100, 100, 500, 400]);
    p1 = loglog(h_values, err_cart, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    hold on;
    p2 = loglog(h_values, err_curv, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    p3 = loglog(h_values, err_curv_scaled, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    
    grid on; box on;
    xlabel('h', 'FontSize', 12);
    ylabel(sprintf('%s Error Norm', error_norm), 'FontSize', 12);
    title(sprintf('%s Error vs Mesh Size (Order %d)', error_norm, order), 'FontSize', 12);
    legend([p1, p2, p3], {'Cartesian', 'Curvilinear (unscaled)', 'Curvilinear (scaled)'}, ...
        'Location', 'best', 'FontSize', 9);
    
%     % Print results
%     fprintf('\n=== ERROR NORM RESULTS (%s) ===\n', error_norm);
%     fprintf('%-15s %12s %12s %12s %12s\n', 'Grid', 'h', 'Cart', 'Curv', 'Curv(sc)');
%     fprintf('%s\n', repmat('-', 55, 1));
%     for g = 1:n_grids
%         fprintf('%-15s %12.4e %12.4e %12.4e %12.4e\n', ...
%             grid_files{sort_idx(g)}, h_values(g), err_cart(g), err_curv(g), err_curv_scaled(g));
%     end
end

function h = compute_mesh_size(GRID, blk, idx)
    i = idx(1); j = idx(2); k = idx(3);
    
    x1 = GRID.gblock(blk).x(i, j, k);
    x2 = GRID.gblock(blk).x(i+1, j, k);
    x3 = GRID.gblock(blk).x(i+1, j+1, k);
    x4 = GRID.gblock(blk).x(i, j+1, k);
    
    y1 = GRID.gblock(blk).y(i, j, k);
    y2 = GRID.gblock(blk).y(i+1, j, k);
    y3 = GRID.gblock(blk).y(i+1, j+1, k);
    y4 = GRID.gblock(blk).y(i, j+1, k);
    
    dx = 0.5 * (sqrt((x2-x1)^2 + (y2-y1)^2) + sqrt((x3-x4)^2 + (y3-y4)^2));
    dy = 0.5 * (sqrt((x4-x1)^2 + (y4-y1)^2) + sqrt((x3-x2)^2 + (y3-y2)^2));
    
    h = max(dx, dy);
end

function plot_convergence_j1_row(grid_files, blk, n_stencil, order, map_order)
    % Plots L2 error convergence for all cells at j=1 across multiple grids
    % Compares Cartesian vs Curvilinear reconstruction
    
    n_grids = length(grid_files);
    
    h_values = zeros(n_grids, 1);
    err_cart = zeros(n_grids, 1);
    err_curv = zeros(n_grids, 1);
    
    parent_dir = pwd;
    
    for g = 1:n_grids
        fprintf('Processing grid %d/%d: %s\n', g, n_grids, grid_files{g});
        
        GRID = load_grid_for_plotting(parent_dir, grid_files{g});
        wall_splines = build_wall_splines(GRID, blk);
        
        block_info_list(1) = block_info_t();
        block_info_list(1).block_id = 1;
        block_info_list(1).Ncells(:) = GRID.gblock(1).Ncells;
        
        imax = GRID.gblock(blk).Ncells(1);
        j_level = 1;
        
        % Use middle cell for h measurement
        i_mid = round(imax / 2);
        h_values(g) = compute_mesh_size(GRID, blk, [i_mid, j_level, 1]);
        
        % Accumulate errors over all i at j=1
        total_error_sq_cart = 0;
        total_error_sq_curv = 0;
        total_volume = 0;
        
        for i = 1:imax
            idx = [i, j_level, 1];
            
            % Build stencil for this cell
            [stencil_idx, ~, ~] = cell_t.build_stencil_alt(blk, idx, n_stencil, block_info_list, true, false);
            
            % Cartesian reconstruction
            [~, xhat_cart] = compute_2D_moments_cartesian(GRID, blk, idx, stencil_idx, order);
            [A_cart, b_cart, col_scale_cart] = assemble_cartesian_ls_system(GRID, stencil_idx, xhat_cart, order);
            coeffs_cart_scaled = A_cart \ b_cart;
            coeffs_cart = coeffs_cart_scaled ./ col_scale_cart';
            
            % Curvilinear reconstruction
            [~, xhat_tn, ~, ~, mapping_data] = compute_2D_moments_curvilinear_poly(GRID, blk, idx, stencil_idx, order, wall_splines, map_order);
            [A_curv, b_curv, col_scale_curv] = assemble_tn_ls_system(GRID, stencil_idx, xhat_tn, order);
            coeffs_curv_scaled = A_curv \ b_curv;
            coeffs_curv = coeffs_curv_scaled ./ col_scale_curv';
            
            % Compute L2 error for this cell only (the reference cell, not whole stencil)
            Q = GRID.gblock(blk).grid_vars.quad(idx(1), idx(2), idx(3));
            x_i = GRID.gblock(blk).grid_vars.cell_c(1, idx(1), idx(2), idx(3));
            y_i = GRID.gblock(blk).grid_vars.cell_c(2, idx(1), idx(2), idx(3));
            
            cell_error_sq_cart = 0;
            cell_error_sq_curv = 0;
            cell_volume = 0;
            
            for q = 1:Q.n_quad
                x_q = Q.quad_pts(1, q);
                y_q = Q.quad_pts(2, q);
                w_q = Q.quad_wts(q);
                
                f_exact = test_function(x_q, y_q);
                
                % Cartesian reconstruction error
                dx = x_q - x_i;
                dy = y_q - y_i;
                f_recon_cart = evaluate_polynomial(coeffs_cart, dx, dy, order);
                cell_error_sq_cart = cell_error_sq_cart + (f_exact - f_recon_cart)^2 * w_q;
                
                % Curvilinear reconstruction error
                [t_q, n_q] = evaluate_curvilinear_mapping(x_q, y_q, mapping_data);
                f_recon_curv = evaluate_polynomial(coeffs_curv, t_q, n_q, order);
                cell_error_sq_curv = cell_error_sq_curv + (f_exact - f_recon_curv)^2 * w_q;
                
                cell_volume = cell_volume + w_q;
            end
            
            total_error_sq_cart = total_error_sq_cart + cell_error_sq_cart;
            total_error_sq_curv = total_error_sq_curv + cell_error_sq_curv;
            total_volume = total_volume + cell_volume;
        end
        
        % Global L2 norm over all j=1 cells
        err_cart(g) = sqrt(total_error_sq_cart / total_volume);
        err_curv(g) = sqrt(total_error_sq_curv / total_volume);
        
        fprintf('  h = %.4e, L2_cart = %.4e, L2_curv = %.4e\n', h_values(g), err_cart(g), err_curv(g));
    end
    
    % Sort by h (descending for proper plotting)
    [h_values, sort_idx] = sort(h_values, 'descend');
    err_cart = err_cart(sort_idx);
    err_curv = err_curv(sort_idx);
    
    % Compute convergence rates
    rate_cart = zeros(n_grids-1, 1);
    rate_curv = zeros(n_grids-1, 1);
    for g = 1:n_grids-1
        rate_cart(g) = log(err_cart(g)/err_cart(g+1)) / log(h_values(g)/h_values(g+1));
        rate_curv(g) = log(err_curv(g)/err_curv(g+1)) / log(h_values(g)/h_values(g+1));
    end
    
    % Plot
    figure('Position', [100, 100, 600, 500]);
    
%     p1 = loglog(h_values, err_cart, 'b-s', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    hold on;
    p1 = loglog(h_values, err_curv, 'r-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % 4th order reference slope
    h_ref = [h_values(1), h_values(end)];
    e_ref_start = err_curv(1) * 2;  % Start slightly above curvilinear
    e_ref = e_ref_start * (h_ref / h_ref(1)).^4;
    p2 = loglog(h_ref, e_ref, 'k--', 'LineWidth', 1.5);
    
    grid on; box on;
    xlabel('h', 'FontSize', 14);
    ylabel('L_2 Error', 'FontSize', 14);
    title(sprintf('Convergence at j=1 (Order %d, Map Order %d)', order, map_order), 'FontSize', 14);
%     legend([p1, p2, p3], {'Cartesian', 'Curvilinear', 'O(h^4) slope'}, 'Location', 'southeast', 'FontSize', 12);
    
    % Print results table
    fprintf('\n=== CONVERGENCE RESULTS (j=1 row) ===\n');
    fprintf('%-12s %12s %12s %8s %12s %8s\n', 'Grid', 'h', 'L2_Cart', 'Rate', 'L2_Curv', 'Rate');
    fprintf('%s\n', repmat('-', 70, 1));
    for g = 1:n_grids
        if g == 1
            fprintf('%-12s %12.4e %12.4e %8s %12.4e %8s\n', ...
                grid_files{sort_idx(g)}, h_values(g), err_cart(g), '-', err_curv(g), '-');
        else
            fprintf('%-12s %12.4e %12.4e %8.2f %12.4e %8.2f\n', ...
                grid_files{sort_idx(g)}, h_values(g), err_cart(g), rate_cart(g-1), err_curv(g), rate_curv(g-1));
        end
    end
    fprintf('\nExpected rate: 4.0\n');
end

function plot_convergence_j11_row(grid_files, blk, n_stencil, order, map_order)
    n_grids = length(grid_files);
    
    h_values = zeros(n_grids, 1);
    err_cart = zeros(n_grids, 1);
    err_curv = zeros(n_grids, 1);
    kappa_cart = zeros(n_grids, 1);
    kappa_curv = zeros(n_grids, 1);
    
    parent_dir = pwd;
    
    for g = 1:n_grids
        fprintf('Processing grid %d/%d: %s\n', g, n_grids, grid_files{g});
        
        GRID = load_grid_for_plotting(parent_dir, grid_files{g});
        wall_splines = build_wall_splines(GRID, blk);
        
        block_info_list(1) = block_info_t();
        block_info_list(1).block_id = 1;
        block_info_list(1).Ncells(:) = GRID.gblock(1).Ncells;
        
        imax = GRID.gblock(blk).Ncells(1);
        j_level = 1;
        
        % Use middle cell for h
        i_mid = round(imax / 2);
        h_values(g) = compute_mesh_size(GRID, blk, [i_mid, j_level, 1]);
        
        total_error_sq_cart = 0;
        total_error_sq_curv = 0;
        total_volume = 0;
        max_kappa_cart = 0;
        max_kappa_curv = 0;
        
        for i = 1:imax
            idx_i = [i, j_level, 1];
            
            % Build stencil for THIS cell
            [stencil_idx_i, ~, ~] = cell_t.build_stencil_alt(blk, idx_i, n_stencil, block_info_list, true, false);
            
            % Cartesian reconstruction
            [~, xhat_cart_i] = compute_2D_moments_cartesian(GRID, blk, idx_i, stencil_idx_i, order);
            [A_cart_i, b_cart_i, col_scale_cart_i] = assemble_cartesian_ls_system(GRID, stencil_idx_i, xhat_cart_i, order);
            coeffs_cart_i = (A_cart_i \ b_cart_i) ./ col_scale_cart_i';
            max_kappa_cart = max(max_kappa_cart, cond(A_cart_i));
            
            % Curvilinear reconstruction
            [~, xhat_tn_i, ~, ~, mapping_data_i] = compute_2D_moments_curvilinear_poly(GRID, blk, idx_i, stencil_idx_i, order, wall_splines, map_order);
            [A_curv_i, b_curv_i, col_scale_curv_i] = assemble_tn_ls_system(GRID, stencil_idx_i, xhat_tn_i, order);
            coeffs_curv_i = (A_curv_i \ b_curv_i) ./ col_scale_curv_i';
            max_kappa_curv = max(max_kappa_curv, cond(A_curv_i));
            
            % Reference point
            x_i = GRID.gblock(blk).grid_vars.cell_c(1, idx_i(1), idx_i(2), idx_i(3));
            y_i = GRID.gblock(blk).grid_vars.cell_c(2, idx_i(1), idx_i(2), idx_i(3));
            
            % Evaluate error at THIS cell only
            Q = GRID.gblock(blk).grid_vars.quad(idx_i(1), idx_i(2), idx_i(3));
            
            cell_error_sq_cart = 0;
            cell_error_sq_curv = 0;
            cell_volume = 0;
            
            for q = 1:Q.n_quad
                x_q = Q.quad_pts(1, q);
                y_q = Q.quad_pts(2, q);
                w_q = Q.quad_wts(q);
                
                f_exact = test_function(x_q, y_q);
                
                % Cartesian
                dx = x_q - x_i;
                dy = y_q - y_i;
                f_recon_cart = evaluate_polynomial(coeffs_cart_i, dx, dy, order);
                
                % Curvilinear
                [t_q, n_q] = evaluate_curvilinear_mapping(x_q, y_q, mapping_data_i);
                f_recon_curv = evaluate_polynomial(coeffs_curv_i, t_q, n_q, order);
                
                cell_error_sq_cart = cell_error_sq_cart + (f_exact - f_recon_cart)^2 * w_q;
                cell_error_sq_curv = cell_error_sq_curv + (f_exact - f_recon_curv)^2 * w_q;
                cell_volume = cell_volume + w_q;
            end
            
            total_error_sq_cart = total_error_sq_cart + cell_error_sq_cart;
            total_error_sq_curv = total_error_sq_curv + cell_error_sq_curv;
            total_volume = total_volume + cell_volume;
        end
        
        err_cart(g) = sqrt(total_error_sq_cart / total_volume);
        err_curv(g) = sqrt(total_error_sq_curv / total_volume);
        kappa_cart(g) = max_kappa_cart;
        kappa_curv(g) = max_kappa_curv;
        
        fprintf('  h=%.2e, κ_cart=%.2e, κ_curv=%.2e, err_cart=%.2e, err_curv=%.2e\n', ...
            h_values(g), kappa_cart(g), kappa_curv(g), err_cart(g), err_curv(g));
    end
    
    % Sort by h descending
    [h_values, sort_idx] = sort(h_values, 'descend');
    err_cart = err_cart(sort_idx);
    err_curv = err_curv(sort_idx);
    kappa_cart = kappa_cart(sort_idx);
    kappa_curv = kappa_curv(sort_idx);
    
    % Compute convergence rates
    fprintf('\n=== CONVERGENCE RESULTS ===\n');
    fprintf('%-12s %10s %10s %8s %10s %8s %10s %10s\n', ...
        'Grid', 'h', 'Cart Err', 'Rate', 'Curv Err', 'Rate', 'κ_cart', 'κ_curv');
    fprintf('%s\n', repmat('-', 90, 1));
    
    for g = 1:n_grids
        if g == 1
            fprintf('%-12s %10.2e %10.2e %8s %10.2e %8s %10.2e %10.2e\n', ...
                grid_files{sort_idx(g)}, h_values(g), err_cart(g), '-', err_curv(g), '-', kappa_cart(g), kappa_curv(g));
        else
            rate_cart = log(err_cart(g-1)/err_cart(g)) / log(h_values(g-1)/h_values(g));
            rate_curv = log(err_curv(g-1)/err_curv(g)) / log(h_values(g-1)/h_values(g));
            fprintf('%-12s %10.2e %10.2e %8.2f %10.2e %8.2f %10.2e %10.2e\n', ...
                grid_files{sort_idx(g)}, h_values(g), err_cart(g), rate_cart, err_curv(g), rate_curv, kappa_cart(g), kappa_curv(g));
        end
    end
    
    % Plot
    figure('Position', [100, 100, 1000, 400]);
    
    % Error plot
    subplot(1,2,1);
    loglog(h_values, err_cart, 'b-s', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    hold on;
    loglog(h_values, err_curv, 'r-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % 4th order slope
    h_ref = [h_values(1), h_values(end)];
    e_ref = err_curv(1) * 2 * (h_ref / h_ref(1)).^4;
    loglog(h_ref, e_ref, 'k--', 'LineWidth', 1.5);
    
    grid on; box on;
    xlabel('h', 'FontSize', 12);
    ylabel('L_2 Error', 'FontSize', 12);
    title('Error Convergence', 'FontSize', 12);
    legend('Cartesian', 'Curvilinear', 'O(h^4)', 'Location', 'southeast');
    
    % Condition number plot
    subplot(1,2,2);
    loglog(h_values, kappa_cart, 'b-s', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    hold on;
    loglog(h_values, kappa_curv, 'r-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    grid on; box on;
    xlabel('h', 'FontSize', 12);
    ylabel('Condition Number', 'FontSize', 12);
    title('Conditioning vs Mesh Size', 'FontSize', 12);
    legend('Cartesian', 'Curvilinear', 'Location', 'northeast');
    
    sgtitle(sprintf('Order %d Reconstruction (j=1 row)', order));
end
function plot_convergence_j1_l1_row(grid_files, blk, n_stencil, order, map_order)
    n_grids = length(grid_files);
    
    h_values = zeros(n_grids, 1);
    err_cart = zeros(n_grids, 1);
    err_curv = zeros(n_grids, 1);
    kappa_cart = zeros(n_grids, 1);
    kappa_curv = zeros(n_grids, 1);
    
    parent_dir = pwd;
    
    for g = 1:n_grids
        fprintf('Processing grid %d/%d: %s\n', g, n_grids, grid_files{g});
        
        GRID = load_grid_for_plotting(parent_dir, grid_files{g});
        wall_splines = build_wall_splines(GRID, blk);
        
        block_info_list(1) = block_info_t();
        block_info_list(1).block_id = 1;
        block_info_list(1).Ncells(:) = GRID.gblock(1).Ncells;
        
        imax = GRID.gblock(blk).Ncells(1);
        j_level = 1;
        
        % Use middle cell for h
        i_mid = round(imax / 2);
        h_values(g) = compute_mesh_size(GRID, blk, [i_mid, j_level, 1]);
        
        total_error_sq_cart = 0;
        total_error_sq_curv = 0;
        total_volume = 0;
        max_kappa_cart = 0;
        max_kappa_curv = 0;
        
        for i = 1:imax
            idx_i = [i, j_level, 1];
            
            % Build stencil for THIS cell
            [stencil_idx_i, ~, ~] = cell_t.build_stencil_alt(blk, idx_i, n_stencil, block_info_list, true, false);
            
            % Cartesian reconstruction
            [~, xhat_cart_i] = compute_2D_moments_cartesian(GRID, blk, idx_i, stencil_idx_i, order);
            [A_cart_i, b_cart_i, col_scale_cart_i] = assemble_cartesian_ls_system(GRID, stencil_idx_i, xhat_cart_i, order);
            coeffs_cart_i = (A_cart_i \ b_cart_i) ./ col_scale_cart_i';
            max_kappa_cart = max(max_kappa_cart, cond(A_cart_i));
            
            % Curvilinear reconstruction
            [~, xhat_tn_i, ~, ~, mapping_data_i] = compute_2D_moments_curvilinear_poly(GRID, blk, idx_i, stencil_idx_i, order, wall_splines, map_order);
            [A_curv_i, b_curv_i, col_scale_curv_i] = assemble_tn_ls_system(GRID, stencil_idx_i, xhat_tn_i, order);
            coeffs_curv_i = (A_curv_i \ b_curv_i) ./ col_scale_curv_i';
            max_kappa_curv = max(max_kappa_curv, cond(A_curv_i));
            
            % Reference point
            x_i = GRID.gblock(blk).grid_vars.cell_c(1, idx_i(1), idx_i(2), idx_i(3));
            y_i = GRID.gblock(blk).grid_vars.cell_c(2, idx_i(1), idx_i(2), idx_i(3));
            
            % Evaluate error at THIS cell only
            Q = GRID.gblock(blk).grid_vars.quad(idx_i(1), idx_i(2), idx_i(3));
            
%             cell_error_sq_cart = 0;
%             cell_error_sq_curv = 0;
%             cell_volume = 0;
%             
%             for q = 1:Q.n_quad
%                 x_q = Q.quad_pts(1, q);
%                 y_q = Q.quad_pts(2, q);
%                 w_q = Q.quad_wts(q);
%                 
%                 f_exact = test_function(x_q, y_q);
%                 
%                 % Cartesian
%                 dx = x_q - x_i;
%                 dy = y_q - y_i;
%                 f_recon_cart = evaluate_polynomial(coeffs_cart_i, dx, dy, order);
%                 
%                 % Curvilinear
%                 [t_q, n_q] = evaluate_curvilinear_mapping(x_q, y_q, mapping_data_i);
%                 f_recon_curv = evaluate_polynomial(coeffs_curv_i, t_q, n_q, order);
%                 
%                 cell_error_sq_cart = cell_error_sq_cart + (f_exact - f_recon_cart)^2 * w_q;
%                 cell_error_sq_curv = cell_error_sq_curv + (f_exact - f_recon_curv)^2 * w_q;
%                 cell_volume = cell_volume + w_q;
%             end

            cell_error_cart = 0;
            cell_error_curv = 0;
            cell_volume = 0;
            total_error_cart = 0;
            total_error_curv = 0;
            n_cells = 0;

            for q = 1:Q.n_quad
                x_q = Q.quad_pts(1, q);
                y_q = Q.quad_pts(2, q);
                w_q = Q.quad_wts(q);

                f_exact = test_function(x_q, y_q);

                % Cartesian
                dx = x_q - x_i;
                dy = y_q - y_i;
                f_recon_cart = evaluate_polynomial(coeffs_cart_i, dx, dy, order);

                % Curvilinear
                [t_q, n_q] = evaluate_curvilinear_mapping(x_q, y_q, mapping_data_i);
                f_recon_curv = evaluate_polynomial(coeffs_curv_i, t_q, n_q, order);

                % L1: absolute value
                cell_error_cart = cell_error_cart + abs(f_exact - f_recon_cart) * w_q;
                cell_error_curv = cell_error_curv + abs(f_exact - f_recon_curv) * w_q;
                cell_volume = cell_volume + w_q;
            end

            % Cell-averaged error (Jalali formula)
            cell_avg_err_cart = cell_error_cart / cell_volume;
            cell_avg_err_curv = cell_error_curv / cell_volume;

            total_error_cart = total_error_cart + cell_avg_err_cart;
            total_error_curv = total_error_curv + cell_avg_err_curv;
            n_cells = n_cells + 1;
            
%             total_error_sq_cart = total_error_sq_cart + cell_error_sq_cart;
%             total_error_sq_curv = total_error_sq_curv + cell_error_sq_curv;
%             total_volume = total_volume + cell_volume;
        end
        
        err_cart(g) = total_error_cart / n_cells;
        err_curv(g) = total_error_curv / n_cells;
        kappa_cart(g) = max_kappa_cart;
        kappa_curv(g) = max_kappa_curv;
        
        fprintf('  h=%.2e, κ_cart=%.2e, κ_curv=%.2e, err_cart=%.2e, err_curv=%.2e\n', ...
            h_values(g), kappa_cart(g), kappa_curv(g), err_cart(g), err_curv(g));
    end
    
    % Sort by h descending
    [h_values, sort_idx] = sort(h_values, 'descend');
    err_cart = err_cart(sort_idx);
    err_curv = err_curv(sort_idx);
    kappa_cart = kappa_cart(sort_idx);
    kappa_curv = kappa_curv(sort_idx);
    
    % Compute convergence rates
    fprintf('\n=== CONVERGENCE RESULTS ===\n');
    fprintf('%-12s %10s %10s %8s %10s %8s %10s %10s\n', ...
        'Grid', 'h', 'Cart Err', 'Rate', 'Curv Err', 'Rate', 'κ_cart', 'κ_curv');
    fprintf('%s\n', repmat('-', 90, 1));
    
    for g = 1:n_grids
        if g == 1
            fprintf('%-12s %10.2e %10.2e %8s %10.2e %8s %10.2e %10.2e\n', ...
                grid_files{sort_idx(g)}, h_values(g), err_cart(g), '-', err_curv(g), '-', kappa_cart(g), kappa_curv(g));
        else
            rate_cart = log(err_cart(g-1)/err_cart(g)) / log(h_values(g-1)/h_values(g));
            rate_curv = log(err_curv(g-1)/err_curv(g)) / log(h_values(g-1)/h_values(g));
            fprintf('%-12s %10.2e %10.2e %8.2f %10.2e %8.2f %10.2e %10.2e\n', ...
                grid_files{sort_idx(g)}, h_values(g), err_cart(g), rate_cart, err_curv(g), rate_curv, kappa_cart(g), kappa_curv(g));
        end
    end
    
    % Plot
    figure('Position', [100, 100, 1000, 400]);
    
    % Error plot
    subplot(1,2,1);
    loglog(h_values, err_cart, 'b-s', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    hold on;
    loglog(h_values, err_curv, 'r-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % 4th order slope
    h_ref = [h_values(1), h_values(end)];
    e_ref = err_curv(1) * 2 * (h_ref / h_ref(1)).^4;
    loglog(h_ref, e_ref, 'k--', 'LineWidth', 1.5);
    
    grid on; box on;
    xlabel('h', 'FontSize', 12);
    ylabel('L_2 Error', 'FontSize', 12);
    title('Error Convergence', 'FontSize', 12);
    legend('Cartesian', 'Curvilinear', 'O(h^4)', 'Location', 'southeast');
    
    % Condition number plot
    subplot(1,2,2);
    loglog(h_values, kappa_cart, 'b-s', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    hold on;
    loglog(h_values, kappa_curv, 'r-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    grid on; box on;
    xlabel('h', 'FontSize', 12);
    ylabel('Condition Number', 'FontSize', 12);
    title('Conditioning vs Mesh Size', 'FontSize', 12);
    legend('Cartesian', 'Curvilinear', 'Location', 'northeast');
    
    sgtitle(sprintf('Order %d Reconstruction (j=1 row)', order));
end