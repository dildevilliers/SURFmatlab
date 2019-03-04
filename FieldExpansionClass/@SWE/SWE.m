classdef SWE < FarFieldExpansion
    
    properties
        Qj
        Qsmn
        NMAX
        MMAX
        r0
        basis
    end
    
    properties (Constant = true, Hidden = true)
        c0 = physconst('Lightspeed');
        eps0 = 8.854187817000001e-12;
        mu0 = 1.256637061435917e-06;
        eta0 = 3.767303134749689e+02;
    end
    
    methods
        
        function obj = SWE(FFobj,NMAXin,MMAXin,r0)
            
            %handle maximum numbers of modes
            %extract farfield vectors from FFobj
            %perform SWE on farfield vectors, get Q-modes and basis
            %Pack each basis mode into its own FarField object
            %Set all properties for object from SWE operation
            
        end
        
        function obj = SWE2FarField(obj,)
            
            %handle maximum numbers of modes
            %extract farfield vectors from FFobj
            %perform SWE on farfield vectors, get Q-modes and basis
            %Pack each basis mode into its own FarField object
            %Set all properties for object from SWE operation
            
        end
        
    end
    
    methods (Static = true)
        obj = FarField2SWE(FFobj,tol,Nbasis,i_FM)
    end
    
end