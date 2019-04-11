classdef (Abstract) FarFieldExpansion
    properties(Abstract, SetAccess = private)
        nBasis
        basis
        nCoeffs
        coeffs
    end
    
    methods
        function plotCoeffs(obj,coeffs)
            %Write me!
        end
        
        function interpolateCoeffs(obj,coeffs,absc,plotflag)
            %Write me!
        end
    end
    
    methods(Abstract, Static)
        coeffs2FarField
        farField2Coeffs
        %getBasisPower
    end
        
end