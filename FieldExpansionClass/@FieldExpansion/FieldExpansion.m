classdef (Abstract) FieldExpansion
    properties 
        nBasis
        nCoeffs
        basis
        coeffs
    end
    
    
    methods (Abstract)
        % FFobj = coeff2field(obj)
        CBFPobj = field2Expansion(FFobj,tol,Nbasis,i_FM);
    end
    
    methods
        % plot coeff
        % plot FF --> to be done by plot
        % plot Basis --> to be done by plot
        % interp coeff
        % select num of coeff
        % select num of basis
    end
end