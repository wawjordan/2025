classdef lsqcon_decomp_t_old
properties
    m (1,1) {mustBeInteger}
    n (1,1) {mustBeInteger}
    p (1,1) {mustBeInteger}
    Q1 (:,:) {mustBeNumeric}
    R1 (:,:) {mustBeNumeric}
    Q2 (:,:) {mustBeNumeric}
    R2 (:,:) {mustBeNumeric}
    Ap (:,:) {mustBeNumeric}
end
methods
    function this = setup(this,A,B)
        % m1 = size(A,1);
        % n1 = size(A,2);
        % m2 = size(B,1);
        this.m = size(A,1);
        this.n = size(A,2);
        this.p = size(B,1);

        % compute the QR factorization B^T=QR
        % DGEQRF
        [this.Q1,Rmat1] = qr(B.');
        
        % solve R(1:m2,1:m2)^T*y = d for y
        % y = R(1:m2,1:m2).' \ d;
        % DTRTRS
        this.R1 = Rmat1(1:this.p,1:this.p).';
        
        % A = A*Q
        % DORMQR
        A_tmp = A*this.Q1;

        % % find z so || A(:,m2+1:n1)*z - ( b - A(:,1:m2)*y ) ||_2 is minimized
        % z = A(:,m2+1:n1) \ ( b - A(:,1:m2)*y );
        % -> QR decomp of A(:,m2+1:n1)
        % DGEQRF
        [this.Q2,this.R2] = qr(A_tmp(:,this.p+1:this.n));
        this.Ap = A_tmp(:,1:this.p);

        % compute y
        %   y = this.R1 \ d;

        % compute z
        %   RHS1 = b - this.Ap * y
        %   get RHS2 from DORMQR with  DTRTRS
        %       RHS2 = this.Q2.' * RHS1
        %       z    = this.R2\RHS2

        % x = Q(:,1:m2)*y + Q(:,m2+1:n1)*z;
        % compute x
        %   x = this.Q2(:,1:p)*y + this.Q2(:,p+1:n)*z;
    end
    function x = solve(this,b,d)
        % compute y
        %   y = this.R1 \ d;

        % compute z
        %   RHS1 = b - this.Ap * y
        %   get RHS2 from DORMQR with  DTRTRS
        %       RHS2 = this.Q2.' * RHS1
        %       z    = this.R2\RHS2

        % x = Q(:,1:m2)*y + Q(:,m2+1:n1)*z;
        % compute x
        %   x = this.Q2(:,1:p)*y + this.Q2(:,p+1:n)*z;
        y = this.R1 \ d;
        RHS1 = b - this.Ap * y;
        RHS2 = this.Q2.' * RHS1;
        z    = this.R2\RHS2;
        x = this.Q1(:,1:this.p)*y + this.Q1(:,this.p+1:this.n)*z;
    end
end
end