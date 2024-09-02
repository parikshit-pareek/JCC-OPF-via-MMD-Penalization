classdef Kernel < handle
    %KERNEL Abstract class for kernels
    
    properties
    end
    
    methods (Abstract)
        % Evaluate this kernel on the data1 and data2 
        Kmat = eval(this, data1, data2);
        
        % Evaluate k(x1, y1), k(x2, y2), .... 
        Kvec = pairEval(this, X, Y);
        
        % Short summary of this kernel. Useful if in the form
        % KerXXX(param1, param2).
        s=shortSummary(this);

        
    end
    
    methods 
        % Just like eval(.) but data1 = data2. This method may be more efficient 
        % for some kernel.
        function Kmat = selfEval(this, data)
            Kmat = this.eval(data, data);
        end

        % Get a random feature map as an object of type FeatureMap from which 
        % random feature approximation of the kernel can be computed. 
        %  - D is a positive integer representing the number of features. 
        %  - d = dimension of the input 
        % Return [] if it is not implemented. 
        %
        function fm = getRandFeatureMap(this, D, d)
            fm = [];
        end
    end

    methods (Static)
        
        %function Ks=candidates(params)
        %    % Create a cell array of Kernel's where params is a cell array
        %    % of parameter p. Each p is what expected from each kernel when
        %    % calling getParam().
        %    % Subclass may override with different number of method arguments. 
        %    error('Subclass of Kernel needs to override this.');
        %end
    end
    
end
