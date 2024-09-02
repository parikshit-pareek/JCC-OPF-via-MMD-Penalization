classdef KGaussian < Kernel
    %KGAUSSIAN Gaussian kernel
    % exp(\| x-y \|_2^2 / (2^sigma2))
    
    properties (SetAccess=private)
        % sigma squared
        sigma2;
    end
    
    methods
        
        function this=KGaussian(sigma2)
            assert(isnumeric(sigma2));
            if length(sigma2) == 1
                this.sigma2 = sigma2;
            else
                for i=1:length(sigma2)
                    this(i) = sigma2(i);
                end
            end
        end
        
        function Kmat = eval(this, X, Y)
            % X, Y are data matrices where each column is one instance
            assert(isnumeric(X));
            assert(isnumeric(Y));

            D2 = bsxfun(@plus, sum(X.^2,1)', sum(Y.^2,1)) - 2*(X'*Y );
            Kmat = exp(-D2./(2*(this.sigma2)));
            
        end
        
        function Kvec = pairEval(this, X, Y)
            assert(isnumeric(X));
            assert(isnumeric(Y));
            
            D2 = sum((X-Y).^2, 1);
            Kvec = exp(-D2./(2*(this.sigma2)));

        end
        
        function s=shortSummary(this)
            s = sprintf('%s(%.3g)', mfilename, this.sigma2);
        end

        function Kmat = selfEval(this, X)
            assert(isnumeric(X));
            s2 = sum(X.^2, 1);

            D2 = bsxfun(@plus, s2, s2') - 2*(X'*X );
            Kmat = exp(-D2./(2*(this.sigma2)));
        end

        % Get a random feature map as an object of type FeatureMap from which 
        % random feature approximation of the kernel can be computed. 
        % D is a positive integer representing the number of features. 
        % Return [] if it is not implemented. 
        %
        function fm = getRandFeatureMap(this, D, d)
            fm = RandFourierGaussMap(this.sigma2, D, d);
        end
    end
    
    methods (Static)

    end
    
    
end

