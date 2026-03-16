classdef monomial_basis_t
    properties
        total_degree (1,1) {mustBeInteger}
        n_dim        (1,1) {mustBeInteger}
        n_terms      (1,1) {mustBeInteger}
        idx          (:,1) {mustBeInteger}
        exponents    (:,:) {mustBeInteger}
        diff_idx     (:,:) {mustBeInteger}
    end
    methods
        function this = monomial_basis_t(n_dim,total_degree)
            if nargin > 0
                this.n_dim        = n_dim;
                this.total_degree = total_degree;
                this.n_terms      = nchoosek( n_dim + total_degree, ...
                                                      total_degree );
                [this.exponents, ...
                 this.idx      , ...
                 this.diff_idx ]  = monomial_basis_t. ...
                                    get_exponents( n_dim, total_degree );
                
            end
        end

        function [val,coef] = eval(this,term,x)
            [val,coef] = monomial_basis_t.evaluate_monomial( ...
                                                       this.n_dim,      ...
                                                       this.exponents, ...
                                                       term, x );
        end

        function [dval,dcoef,coef] = deval(this,term,x,order)
            [dval,dcoef,coef] = monomial_basis_t.                       ...
                                evaluate_monomial_derivative(...
                                                       this.n_dim,      ...
                                                       this.exponents,  ...
                                                       term, x, order );
        end

        function check_gradient_idexing(this)
            tuple_str = @(list) sprintf([repmat('%d,',1,numel(list)-1),...
                                         '%d'],list);
            for d = this.total_degree-1:-1:1
                for term = this.idx(d)+1:this.idx(d+1)
                    grad_idx = this.diff_idx(:,term);
                    str = [ tuple_str([d,term]), ' : [',                ...
                            tuple_str(this.exponents(:,term)), '] : (', ...
                            tuple_str(grad_idx), ')' ];
                    fprintf('%s\n',str);
                end
            end
            term = 1;
            grad_idx = this.diff_idx(:,term);
            str = [ tuple_str([d,term]), ' : [',                        ...
                    tuple_str(this.exponents(:,term)), '] : (',         ...
                    tuple_str(grad_idx), ')' ];
            fprintf('%s\n',str);
        end

    end
    methods (Static)
        function [val,coef] = evaluate_monomial( n_dim,     ...
                                                 exponents, ...
                                                 term,      ...
                                                 x )
            val = 1.0;
            coef = 1;
            for d = 1:n_dim
                for i = exponents(d,term):-1:1
                    val = val*x(d);
                    coef = coef * i;
                end
            end
        end

        function [dval,dcoef,coef] = evaluate_monomial_derivative(      ...
                                                             n_dim,     ...
                                                             exponents, ...
                                                             term,      ...
                                                             x,         ...
                                                             order )
            dval  = 0.0;
            dcoef = 1;
            coef  = 1;
            if ( any( exponents(:,term) - order(:) < 0 ) ); return; end
            dval = 1.0;
            for d = 1:n_dim
                for i = exponents(d,term):-1:exponents(d,term)-order(d)+1
                    dcoef = dcoef * i;
                end
                for i = exponents(d,term)-order(d):-1:1
                    coef = coef * i;
                    dval = dval*x(d);
                end
            end
        end

        function [exponents,idx,diff_idx] = get_exponents( n_dim, degree )
            n_terms   = nchoosek( n_dim + degree, degree );
            exponents = zeros(n_dim,n_terms);
            idx       = zeros(degree+1,1);
            diff_idx  = zeros(n_dim,n_terms);
            cnt = 0;
            for curr_total_degree = 0:degree
                nsub(1:n_dim) = curr_total_degree + 1;
                N_full_terms = (curr_total_degree+1)^n_dim;
                for j = 0:N_full_terms
                    tmp_exp = monomial_basis_t.global2local(j+1,nsub)-1;
                    if ( sum(tmp_exp) == curr_total_degree )
                        cnt = cnt + 1;
                        exponents(:,cnt) = tmp_exp;
                    end
                end
                idx(curr_total_degree+1) = cnt;
            end
            diff_idx(1:n_dim,1:n_terms) = -1;
            if ( degree == 0 ); return; end
            for j = 1:idx(degree)
                tmp_exp = exponents(:,j);
                curr_total_degree = sum(tmp_exp);
                cnt = 0;
                for i = idx(curr_total_degree+1)+1:idx(curr_total_degree+2)
                    if ( sum( abs(exponents(:,i) - tmp_exp) ) == 1 )
                        cnt = cnt +1;
                        diff_idx(cnt,j) = i;
                    end
                end
            end
        end

        function iSub = global2local(iG,nSub)
            nDims = numel(nSub);
            iSub = zeros(1,nDims);
            if (nDims==1)
                iSub(1) = iG;
                return
            end
            p = prod(nSub);
            iGtmp = iG;
            for i = nDims:-1:1
                p = fix( p/nSub(i) );
                iTmp = mod(iGtmp-1,p)+1;
                iSub(i) = fix( (iGtmp-iTmp)/p ) + 1;
                iGtmp = iTmp;
            end
        end

        function iG = local2global(iSub,nSub)
            iSub = iSub(:);
            nSub = nSub(:);
            nDims = numel(iSub);
            p = 1;
            iG = 1;
            for i = 1:nDims
                iG = iG + ( iSub(i) - 1 )*p;
                p = p*nSub(i);
            end
        end
    end

end
