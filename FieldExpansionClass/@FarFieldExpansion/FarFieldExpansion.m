classdef (Abstract) FarFieldExpansion
    properties(Abstract)
        nBasis
        basis
        nCoeffs
        coeffs
        nCoeffs
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
        expansion2FarField
        farField2Coefficients
        %getBasisPower
    end
    
%     methods
%         function plotCoeffs(obj,coeffs)
%             %Write me!
%         end
%         
%         function interpolateCoeffs(obj,coeffs,absc,plotflag)
%             %Write me!
%         end
%     end
        
end