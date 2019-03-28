classdef (Abstract) FarFieldExpansion
    properties(Abstract)
        nBasis
        basis
        coeffs
    end
    
    methods
        function plotCoeffs(obj,coeffs)
            %Write me!
        end
        
        function interpolateCoeffs(obj,coeffs,absc)
            %Write me!
        end
    end
    
    methods(Abstract, Static)
        expansion2FarField
        farField2Expansion
    end
    
end