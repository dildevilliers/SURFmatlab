classdef (Abstract) FarFieldExpansion
    properties(Abstract)
        nBasis
        basis
    end
    
    methods(Abstract, Static)
        expansion2FarField
        farField2Expansion
    end
    
end