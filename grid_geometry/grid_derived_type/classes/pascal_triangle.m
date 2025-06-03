classdef pascal_triangle
    properties
        lut (:,1)  {mustBeInteger}
        % lut = [ 1                  ], // n=0
        %       [ 1, 1               ], // n=1
        %       [ 1, 2, 1            ], // n=2
        %       [ 1, 3, 3, 1         ], // n=3
        %       [ 1, 4, 6, 4, 1      ], // n=4
        %       ...
        Nmax (1,1) {mustBeInteger} = 11
    end

    methods
        function this = pascal_triangle(Nmax)
            if nargin > 0
                this.Nmax = Nmax+1;
            end
            this.lut = this.generate_pascal_triangle_lut(this.Nmax);
        end
        function val = binomial_coefficient(this,n,k)
            if k > n+1
                error('K must be an integer between 0 and N.')
            end
            if (n+1>this.Nmax)||(k+1>this.Nmax)
                error('Error, insufficient Nmax in LUT.')
            end
            val = this.index_packed_ut(this.lut,k+1,n+1);
        end
    end

    methods (Static)
        function p = generate_pascal_triangle_lut(n)
            p = ones( ( (n+1)*(n+2) )/2,1);
            for j = 3:n+1
                for i = 2:j-1
                    k0 =   i   + ( (   j  *(j-1) )/2 );
                    k1 =   i   + ( ( (j-1)*(j-2) )/2 );
                    k2 = (i-1) + ( ( (j-1)*(j-2) )/2 );
                    p(k0) = p(k1) + p(k2);
                end
            end
        end
        function val = index_packed_ut(p,i,j)
            val = p(i+(j*(j-1))/2);
        end
    end
end