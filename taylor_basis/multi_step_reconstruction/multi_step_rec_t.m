classdef multi_step_rec_t

    properties
        n_nbors  (1,1) {mustBeInteger} = 1
        idx      (:,1) {mustBeInteger}
        nbor_idx (:,:) {mustBeInteger}
        nbor_face_id (:,1) {mustBeInteger}
        basis    (1,1) zero_mean_basis2
        quad     (1,1) quad_t
        fquad    (:,1) quad_t
        fvec     (:,1) face_vec_t
        n_vars    (1,1) {mustBeInteger}
        avgs     (:,1)   {mustBeNumeric}
        coefs    (:,:)   {mustBeNumeric}
        rec_LHS_inv  (:,1) cell
        rec_LHS_tmp  (:,1) cell
        rec_LHS      (:,1) cell
        rec_RHS  (:,1) cell
    end

    methods
        function this = multi_step_rec_t(nbor_idxs,nbor_face_id,T,Q,FQ,FN,n_vars)
            if nargin<1
                return
            end
            this.idx      = nbor_idxs(:,1);
            this.nbor_idx = nbor_idxs(:,2:end);
            this.nbor_face_id = nbor_face_id;
            this.n_nbors  = size(nbor_idxs,2) - 1;
            this.basis    = T;
            this.quad     = Q;
            this.fquad    = [FQ{:}];
            this.fvec     = [FN{:}];
            this.n_vars   = n_vars;
            this.avgs     = zeros(n_vars,1);
            this.coefs    = zeros(this.basis.n_terms,n_vars);
        end

        function val = eval_reconstruction(this,x,n_terms)
            val = zeros(this.n_vars,1);
            Bx = zeros(this.basis.n_terms,1);
            if nargin > 2
                for i = 1:n_terms
                    Bx(i) = this.basis.eval(i,x);
                end
            else
                for i = 1:this.basis.n_terms
                    Bx(i) = this.basis.eval(i,x);
                end
            end

            for v = 1:this.n_vars
                val(v) = sum(Bx.*this.coefs(:,v));
            end
        end

        function this = set_cell_avg(this,fun)
            n_dim  = this.basis.n_dim;
            n_quad = this.quad.n_quad;
            volume = this.quad.integrate(ones(1,n_quad));
            tmp_val = zeros(1,n_quad);
            for v = 1:numel(fun)
                for n = 1:n_quad
                   tmp_x      = num2cell(this.quad.quad_pts(1:n_dim,n));
                   tmp_val(n) = fun{v}(tmp_x{:});
                end
                this.avgs(v) = this.quad.integrate(tmp_val)/volume;
            end
            this.coefs(1,:) = this.avgs;
        end

        function this = setup_rec(this,nbors)
            % this.self_LHS = this.get_self_LHS(n1,nbors);
        end

        function this = get_LHS_tmp(this,nbors,step)
            col_start = 1;
            col_end   = zero_mean_basis2.get_n_terms(this.basis.n_dim+step,step) - 1;

            if (step == 1)
                BM = this.multi_step_rec_t.get_basis_matrix(nbors);
                this.rec_LHS_inv{step} = pinv( BM(:,col_start:col_end) );
                this.rec_LHS_tmp{step} = this.rec_LHS_inv{step} * BM;
            else
                this.rec_LHS_inv{step} = pinv( this.rec_LHS{step-1}(:,col_start:col_end) );
                this.rec_LHS_tmp{step} = this.rec_LHS_inv{step} * this.rec_LHS{step-1};
            end
            % A+_{i,1} * A_i
        end

        function this = assemble_LHS(this,nbors,step)
            row_sz = size(this.rec_LHS_inv{step},1);
            col_sz = this.basis.n_terms - 1;
            for n = 1:numel(nbors)
                row_sz = row_sz + size(nbors(n).rec_LHS_inv{step},1);
            end
            this.rec_LHS{step} = zeros(row_sz,col_sz);

            row_start = 1;
            row_end   = size(this.rec_LHS_inv{step},1);
            this.rec_LHS{step}(row_start:row_end,:) = this.rec_LHS_tmp{step};
            % A+_{i,1} * A_i
            for n = 1:numel(nbors)
                row_start = row_end + 1;
                row_end   = row_end + size(nbors(n).rec_LHS_inv{step},1);
                T = multi_step_rec_t.get_t_matrix(this,nbors(n));
                this.rec_LHS{step}(row_start:row_end,:) = this.rec_LHS_tmp{step} * T;
                % A+_{jn,1} * A_jn
            end
        end

        % function this = get_next_LHS(this,nbors,step)
        %     if (step == 1)
        %         BM = multi_step_rec_t.get_basis_matrix(this,nbors);
        %         this = this.get_LHS_inv1(nbors);
        %         this.rec_LHS{1} = this.rec_LHS_inv{1} * BM;
        % 
        %         sz1 = size(A1);
        %     end
        % end
    end
    methods (Static)

        function CELLS = setup_reconstructions(CELLS,steps)
            for i = 1:steps
                for n = 1:numel(CELLS)
                    CELLS(n) = CELLS(n).
        end

        function BM = get_basis_matrix(CELL,NBORS)
            % evaluate basis functions
            BM = zeros(CELL.n_nbors,CELL.basis.n_terms);
            for j = 1:CELL.basis.n_nbors
                for i = 2:CELL.basis.n_terms
                    BM(i-1,j) = NBORS(j).basis.eval(i,CELL.basis.x_ref);
                end
            end
        end

        function T = get_t_matrix(CELL,NBOR)
            T = zeros(CELL.basis.n_terms-1,CELLS.basis.n_terms-1);

            % row & column scaling
            sj  = zeros(CELL.basis.n_terms-1,1);
            si  = zeros(CELL.basis.n_terms-1,1);
            for j = 2:CELL.basis.n_terms
                [h,f] = CELL.basis.terms(j).eval(NBOR.basis.h_ref);
                sj(j-1) = f/h;
                [h,f] = CELL.basis.terms(j).eval(NBOR.basis.h_ref);
                si(j-1) = f/h;
            end

            dij = NBOR.basis.x_ref - CELL.basis.x_ref;
            
            for j = 2:CELL.basis.n_terms
                for i = 2:j
                    order = CELL.basis.terms(i).exponents;
                    dval = CELL.basis.terms(j).deval( dij, order );
                    T(i-1,j-1) = (1/sj(i-1)) * dval * si(j-1);
                end
            end
        end

        % function LHS = get_sub_matrix(CELL,NBORS,iter)
        % 
        %     if iter == 0
        %         BM = multi_step_rec_t.get_basis_matrix(CELL,NBORS,2,nchoosek(CELL.basis.n_dim + 1,1));
        %         LHS = 
        %     else
        %         LHS = get_sub_matrix(CELL,NBORS,iter-1);
        % 
        % 
        % end

        % function CELLS = setup_all(n1,CELLS)
        %     n_dim   = CELLS(1).basis.n_dim;
        %     dim_sz  = num2cell(1:n_dim);
        %     sz      = size(CELLS,dim_sz{:});
        %     n_cells = numel(CELLS);
        %     for i = 1:n_cells
        %         n_nbor = CELLS(i).n_nbors;
        %         nbor_id = zeros(n_nbor,1);
        %         for n = 1:n_nbor
        %             nbor_id(n) = zero_mean_basis.local_to_global(CELLS(i).nbor_idx(1:n_dim,n),sz);
        %         end
        %         CELLS(i) = CELLS(i).setup_rec(n1,CELLS(nbor_id));
        %     end
        % end

        % function [CELLS,LHS,RHS,coefs] = perform_reconstruction(n1,CELLS)
        %     n_terms_local = CELLS(1).basis.n_terms - n1;
        %     n_vars        = CELLS(1).n_vars;
        %     n_cells       = numel(CELLS);
        %     LHS = multi_step_rec_t.construct_full_LHS_matrix(n1,CELLS);
        %     RHS = multi_step_rec_t.construct_full_RHS(n1,CELLS);
        %     coefs = zeros(n_terms_local*n_cells,CELLS(1).n_vars);
        %     for v = 1:n_vars
        %         coefs(:,v) = LHS\RHS(:,v);
        %     end
        % 
        %     for i = 1:n_cells
        %         row_start = (i-1)*n_terms_local + 1;
        %         row_end   =   i  *n_terms_local;
        %         n_terms   = CELLS(i).basis.n_terms;
        %         % CELLS(i).coefs(1,:) = CELLS(i).avgs;
        %         CELLS(i).coefs(n1+1:n_terms,:) = coefs(row_start:row_end,:);
        %     end
        % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [X,Fi,Fj,Qi] = evaluate_taylor_basis_and_derivatives(CELL,NBORS,n1,npts)

            n_nbor  = CELL.n_nbors;
            n_dim   = CELL.basis.n_dim;
            n_deriv = CELL.basis.n_terms;
            n_terms = CELL.basis.n_terms - n1;

            % cell array for evaluation grids
            X = cell(n_nbor+1,1);

            X{1} = get_local_interp_grid(CELL.quad,n_dim,npts);
            for k = 1:n_nbor
                X{k+1} = get_local_interp_grid(NBORS(k).quad,n_dim,npts);
            end

            % cell arrays for evaluated basis functions + derivatives
            Fi = cell(n_nbor,1);
            Fj = cell(n_nbor,1);

            Qi = cell(n_nbor,1);
            % cell
            for k = 1:n_nbor
                % basis functions
                Fi{k} = cell(n_terms,1);
                Fj{k} = cell(n_terms,1);
                for j = 1:n_terms
                    % derivatives
                    Fi{k}{j} = cell(n_deriv,1);
                    Fj{k}{j} = cell(n_deriv,1);
                end

                % face quad_pts
                Qi{k} = CELL.fquad( CELL.nbor_face_id(k) ).quad_pts(1:2,:).';
            end

            for i = 1:n_deriv
                for j = 1:n_terms
                    for k = 1:n_nbor
                        dij = CELL.basis.x_ref - NBORS(k).basis.x_ref;
                        dij_mag = sqrt( sum(dij.^2));
                        Fi{k}{j}{i} = arrayfun( @(x1,y1)    CELL.basis.calc_basis_derivative(j,i,[x1;y1],dij_mag), X{1}{1}, X{1}{2} );
                        Fj{k}{j}{i} = arrayfun( @(x1,y1)NBORS(k).basis.calc_basis_derivative(j,i,[x1;y1],dij_mag), X{k+1}{1}, X{k+1}{2} );
                    end
                end
            end

        end

        function X = get_local_interp_grid(Q,n_dim,npts)
            % get number of quadrature points (per dimension)
            n_quad = repmat( round( Q.n_quad.^(1/n_dim) ), 1, n_dim );
            quad_ref = quad_t.create_quad_ref_2D( n_quad(1) );


            n_plot_pts = ones(3,1);
            n_plot_pts(1:n_dim) = npts;

            xtmp = reshape(Q.quad_pts(1,:),n_quad);
            ytmp = reshape(Q.quad_pts(2,:),n_quad);
            ztmp = reshape(Q.quad_pts(2,:),n_quad);

            L = lagrange_interpolant(n_quad(1) + 1);
            % extrapolate the parametric space from the innermost quadrature points
            % Note: needs to be a consistent order of quadrature for the geometry representation
            [x1,y1,z1] = L.map_grid_from_quad(xtmp,ytmp,ztmp,n_plot_pts(1),...
                n_plot_pts(2),...
                n_plot_pts(3), quad_ref );
            X = {x1,y1,z1};
        end

        function macro_plot_self(fig_h,X,Fi,nbor,term,deriv,varargin)
            set(0, 'currentfigure', fig_h);
            hold on;
            surf(X{1}{1},X{1}{2},Fi{nbor}{term}{deriv},varargin{:})
            view(-14,12)
        end

        function macro_plot(fig_h,X,Fj,Qi,nbor,term,deriv,varargin1,varargin2)
            set(0, 'currentfigure', fig_h);
            hold on;
            surf(X{1+nbor}{1},X{1+nbor}{2},Fj{nbor}{term}{deriv},varargin1{:})
            % plot where the quadrature points are
            plot(Qi{nbor}(:,1),Qi{nbor}(:,2),varargin2{:})
            view(-14,12)
        end

        % function macro_plot(fig_h,X,Fi,Fj,Qi,nbor,term,deriv)
        %     set(0, 'currentfigure', fig_h);
        %     hold on;
        %     % plot where the quadrature points are
        %     plot(Qi{nbor}(:,1),Qi{nbor}(:,2),'ko')
        %     surf(X{1}{1},X{1}{2},Fi{nbor}{term}{deriv},'FaceColor','b','EdgeColor','k')
        %     surf(X{1+nbor}{1},X{1+nbor}{2},Fj{nbor}{term}{deriv},'FaceColor','r','EdgeColor','k')
        %     view(-14,12)
        % end

        function macro_plot_quads(fig_h,Qi,nbor,varargin)
            set(0, 'currentfigure', fig_h);
            % plot where the quadrature points are
            plot(Qi{nbor}(:,1),Qi{nbor}(:,2),varargin{:})
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end