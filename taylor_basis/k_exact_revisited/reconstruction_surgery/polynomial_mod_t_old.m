classdef polynomial_mod_t_old
    properties
        degree     (1,1) {mustBeInteger} = 0    % total degree of polynomial
        n_dim      (1,1) {mustBeInteger} = 0    % # of spatial dimensions
        n_var      (1,1) {mustBeInteger} = 0    % # of variables
        n_terms    (1,1) {mustBeInteger} = 0    % # of terms in the polynomial
        exponents  (:,:) {mustBeInteger}        % exponents
        coeff      (:,:) {mustBeNumeric}        % coefficients
        x_ref      (:,1) {mustBeNumeric}        % reference location
    end
    methods
        function this = read_from_file(this,filename)
            % open the file
            fid = fopen(filename,"r");

            % read # of vars and # of terms
            tmp = fscanf(fid,'%d',3);
            this.n_dim   = tmp(1);
            this.n_var   = tmp(2);
            this.n_terms = tmp(3);

            this.x_ref = zeros(this.n_dim,1);
            this.exponents = zeros(this.n_dim,this.n_terms);
            this.coeff     = zeros(this.n_var,this.n_terms);

            % read x_ref
            x_ref_tmp = fscanf(fid,'%f',this.n_dim);
            this.x_ref = zeros(3,1);
            this.x_ref(1:this.n_dim) = x_ref_tmp(:);

            % read exponents
            for i = 1:this.n_terms
                this.exponents(:,i) = fscanf(fid,'%d',this.n_dim);
            end

            this.degree = max(this.exponents,[],"all");

            % read coefficients
            for j = 1:this.n_var
                for i = 1:this.n_terms
                    this.coeff(j,i) = fscanf(fid,'%f',1);
                end
            end

            % close the file
            fclose(fid);
        end
        function this = polynomial_mod_t_old( n_var, n_dim, degree )
            if nargin ~= 3
                return
            else
                this.degree    = degree;
                this.n_dim     = n_dim;
                this.n_var     = n_var;
                this.n_terms   = nchoosek(n_dim+degree,degree);
                % this.exponents = zeros(n_dim,this.n_terms);
                this.coeff     = zeros(n_var,this.n_terms);
                this.x_ref     = zeros(n_dim,1);
                this.exponents = pexponents_alt_old( n_dim, degree );
            end
        end

        function poly_eval = evaluate( this, x_eval )
            poly_eval = zeros(this.n_var,1);
            x = x_eval(:) - this.x_ref;
            for n = 1:this.n_terms
                % Evaluate the Corresponding Monomial
                % monomial_eval = 1.0;
                % for d = 1:this.n_dim
                %     monomial_eval = monomial_eval * ( x(d)^this.exponents(d,n) );
                % end
                monomial_eval = evaluate_monomial_d0_old(x,this.exponents(:,n),this.n_dim);
                % Add Contribution to Polynomial Evaluation
                poly_eval = poly_eval + this.coeff(:,n)*monomial_eval;
            end
        end

        function poly_eval = evaluate_limited( this, x_eval, mask )
            poly_eval = zeros(this.n_var,1);
            x = x_eval(:) - this.x_ref;
            for n = 1:this.n_terms
                % Evaluate the Corresponding Monomial
                % monomial_eval = 1.0;
                % for d = 1:this.n_dim
                %     monomial_eval = monomial_eval * ( x(d)^this.exponents(d,n) );
                % end
                monomial_eval = evaluate_monomial_d0_old(x,this.exponents(:,n),this.n_dim);
                % Add Contribution to Polynomial Evaluation
                poly_eval = poly_eval + mask(n)*this.coeff(:,n)*monomial_eval;
            end
        end

        function grad_eval = gradient( this, x_eval )
            grad_eval = zeros(this.n_var,this.n_dim);
            x = x_eval - this.x_ref;
            for n = 1:this.n_terms
                % % partial derivatives in dimension d
                % for d = 1:this.n_dim
                %     d_exponents = this.exponents(:,n);
                %     d_exponents(d) = d_exponents(d)-1;
                %     monomial_eval = 1.0;
                %     for di = 1:this.n_dim
                %         for i = d_exponents(di):-1:1
                %             monomial_eval = monomial_eval*x(di);
                %         end
                %     end
                %     grad_eval(:,d) = grad_eval(:,d) + this.exponents(d,n)*this.coeff(:,n)*monomial_eval;
                % end
                monomial_eval = evaluate_monomial_d1_old(x,this.exponents(:,n),this.n_dim);
                for d = 1:this.n_dim
                    grad_eval(:,d) = grad_eval(:,d) + this.coeff(:,n)*monomial_eval(d);
                end
            end
        end
        function hess_eval = hessian( this, x_eval )
            sz = nchoosek(this.n_dim+1,2);
            hess_eval = zeros(this.n_var,sz);
            x = x_eval - this.x_ref;
            for n = 1:this.n_terms
                monomial_eval = evaluate_monomial_d2_old(x,this.exponents(:,n),this.n_dim);
                for d = 1:sz
                    hess_eval(:,d) = hess_eval(:,d) + this.coeff(:,n)*monomial_eval(d);
                end
            end
        end
        function grad_eval = gradient_old( this, x_eval )
            grad_eval = zeros(this.n_var,this.n_dim);
            x = x_eval - this.x_ref;
            % Compute Gradients
            for vi = 1:this.n_var
                % Compute Gradient for a Given Variable
                var_grad = zeros(this.n_dim,1);
                for n = 1:this.n_terms
                    % Extract Term Coefficient & Exponent
                    term_coeff    = this.coeff(vi,n);
                    term_exponent = this.exponents(:,n);
                    for di = 1:this.n_dim
                        % Determine if Current Term is a Function of the Current Dimension
                        if ( term_exponent(di) == 0 )
                            continue;
                        end
                        % Derivative Coefficient
                        d_coeff = term_coeff*term_exponent(di);
                        % Derivative Monomial Exponent
                        d_exponent     = term_exponent;
                        d_exponent(di) = d_exponent(di) - 1;
                        % Evaluate Monomial
                        monomial_eval = 1.0;
                        for d = 1:this.n_dim
                            monomial_eval = monomial_eval * ( x(d)^d_exponent(d) );
                        end
                        % monomial_eval = evaluate_monomial( x_eval, this.x_ref, ...
                        %                                    d_exponent, this.n_dim==2 );
                        % Add Contribution
                        var_grad(di) = var_grad(di) + d_coeff*monomial_eval;
                    end
                end
                grad_eval(vi,:) = var_grad;
            end
        end
        function linear_coeffs = compute_linearization(this,x_eval)
            linear_coeffs = zeros(this.n_terms,1);
            delta = x_eval(1:this.n_dim) - this.x_ref;
            for n = 1:this.n_terms
                linear_coeffs(n) = evaluate_monomial_d0_old( delta, this.exponents(:,n), this.n_dim );
            end
        end
    end

end