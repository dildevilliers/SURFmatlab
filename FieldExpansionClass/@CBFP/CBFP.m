classdef CBFP < FieldExpansion 
    properties
        UR
        S
        VR
        sigma_n
    end
    
        methods
            function obj = CBFP(nBasis,nCoeffs,UR,S,VR,basis,coeffs,sigma_n)
                % function obj = CBFP(nBasis,S,U,V,R,basisFuncs,coeffs,sigma_th_n,sigma_ph_n)
                % CBFP object constructor method
                %
                % Inputs:   nBasis  - Number of basis functions
                %           nCoeffs - Number of coefficients
                %           UR      - SVD: reduced left singular vectors
                %           S      - Singular values
                %           VR      - SVD: reduced right singular vectors
                %           basis   - Set of orthonormal basis functions
                %           coeffs  - Coefficients
                %           sigma_n - Normalised singular values {diag(S)./S(1,1)}
                %
                % Output:   obj     - CBFP object
                
                obj.nBasis = nBasis;
                obj.nCoeffs = nCoeffs;
                obj.UR = UR;
                obj.S = S;
                obj.VR = VR;
                obj.basis = basis;
                obj.coeffs = coeffs;
                obj.sigma_n = sigma_n;
            end
        end
    
        methods (Static = true)
            CBFPobj = field2Expansion(FFobj,tol,Nbasis,i_FM)
        end
    
end