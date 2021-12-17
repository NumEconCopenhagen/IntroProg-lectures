classdef classier
    %CLASSIER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        A = pi; 
    end
    
    methods (Static)
        function f = f1(x,y)
            f = max(exp(x),classier.A*y); 
        end; 
        
        function f = f2(x,y,z)
            f = x + 3*x*classier.f1(y,z); 
        end; 
    end
    
end

