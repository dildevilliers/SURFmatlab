classdef CBFP < FarFieldExpansion 
    properties
        tol
        i_FM 
        UR
        SR
        VR
        sigma_n
    end
    
    properties (SetAccess = immutable, Hidden = true)
        isFreqBasis
        freqRange
    end
    
        methods
            function obj = CBFP(FFobj,tol,i_FM,nBasis)
                % function obj = CBFP(FFobj,tol,i_FM,nBasis)
                % CBFP object constructor method
                %
                % Inputs:   FFobj   - Farfield object
                %           tol     - Tolerance 
                %           i_FM    - Struct of indices for CBFP support points
                %                     i_FM.f => indices across frequency
                %                     i_FM.x => indices across the design space
                %           nBasis  - Number of basis functions
                %
                % Output:   CBFPobj - CBFP object
                
                if nargin == 0
                    % Expand gausian beam pattern
                else
                    % INPUT CHECKS
                    
                    if nargin < 2
                        tol = 1e-100;
                        i_FM = [];
                        nBasis = inf;
                    elseif nargin < 3 
                        i_FM = [];
                        nBasis = inf;
                    elseif nargin < 4 
                        nBasis = inf;
                    end

                    % Ignore some function inputs
                    if isempty(tol), tol = 1e-100; end
                    if isempty(i_FM), i_FM = []; end
                    if isempty(nBasis), nBasis = inf; end
                    
                    nFFobjs = length(FFobj);
                    
                    if nFFobjs == 1
                        
                        isFreqBasis = 1;
                        
                        if isfield(i_FM,'x'), warning('No geometric variation. i_FM.x ignored'), end
                        if isfield(i_FM,'f'),assert(length(i_FM.f) <= FFobj.Nf,'Specified number of indices,i_FM.f larger than number of frequency samples,FFobj.Nf'), end
                        if isfield(i_FM,'f'),assert(max(i_FM.f) <= FFobj.Nf,'Largest specified index in i_FM.f larger than maximum possible index FFobj.Nf'), end
                        
                        if nBasis ~= inf
                            if isfield(i_FM,'f'), assert(nBasis <= FFobj.Nf,'nBasis should not exceed FFobj.Nf'), end
                            if isfield(i_FM,'f'), assert(nBasis <= length(i_FM.f),'nBasis should not exceed the total number of possible basis functions, as selected by i_FM.f'), end
                        end
                    else
                        Nf = FFobj(1,1).Nf;
                        if Nf == 1
                            if isfield(i_FM,'f'), warning('Frequency treated as a normal parameter. i_FM.f ignored'), end
                            if isfield(i_FM,'x'),assert(length(i_FM.x) <= nFFobjs,'Specified number of indices,i_FM.x larger than number of geometric variations, num of FFobjs'), end
                            if isfield(i_FM,'x'),assert(max(i_FM.x) <= nFFobjs,'Largest specified index in i_FM.x larger than maximum possible index, num of FFobjs'), end
                            
                            if nBasis ~= inf
                                assert(nBasis <= nFFobjs,'nBasis should not exceed the maximum possible number of basis functions')
                                if isfield(i_FM,'x'), assert(nBasis <= length(i_FM.x),'nBasis should not exceed the total number of possible basis functions, as selected by i_FM.x'), end
                            end
                        else
                            if isfield(i_FM,'f'),assert(length(i_FM.f) <= Nf,'Specified number of indices,i_FM.f larger than number of frequency samples,FFobj.Nf'), end
                            
                            if isfield(i_FM,'f'),assert(max(i_FM.f) <= Nf,'Largest specified index in i_FM.f larger than maximum possible index FFobj.Nf'), end
                            
                            if isfield(i_FM,'x'),assert(length(i_FM.x) <= nFFobjs,'Specified number of indices,i_FM.x larger than number of geometric variations, num of FFobjs'), end
                            if isfield(i_FM,'x'),assert(max(i_FM.x) <= nFFobjs,'Largest specified index in i_FM.x larger than maximum possible index, num of FFobjs'), end
                            
                            if nBasis ~= inf
                                assert(nBasis <= nFFobjs,'nBasis should not exceed the maximum possible number of basis functions')
                                if isfield(i_FM,'x'), assert(nBasis <= length(i_FM.x),'nBasis should not exceed the total number of possible basis functions, as selected by i_FM.x'), end
                            end
                        end
                    end
                                  
                    %% FREQUENCY ONLY VARIATION
                    nFFobjs = length(FFobj);
                    if nFFobjs == 1
                        % Extract useful parameters
                        Nf = FFobj.Nf;
                        E1 = FFobj.E1;
                        E2 = FFobj.E2;
                        Nang = FFobj.Nang;
                        
                        % Determine FF indices from which CBFP matrix is built                       
                        if isfield(i_FM,'f')
                            range_f = reshape(i_FM.f,1,length(i_FM.f));
                        else
                            range_f = 1:Nf;
                        end
                        
                        % Get the CBFP matrix
                        FM_E1 = [];
                        FM_E2 = [];
                        FM_E1 = [FM_E1 E1(:,range_f)];
                        FM_E2 = [FM_E2 E2(:,range_f)];
                        
                        FM = [FM_E1;FM_E2];
    
                        % Get the SVD of the FM matrix
                        [U,S,V] = svd(FM,'econ');
                        
                        Snorm = S(1,1);
                        sigma_n = diag(S)./Snorm;   % Normalization factor for sigma - this leaves an option to zero one of the field components if required and get rid of noise

                        % Sort out the SVD matrix sizes due to reduced number of significant singular values
                        if nBasis < inf
                            NR = nBasis;
                        else
                            NR = length(find(sigma_n > tol));
                        end
                        
                        if NR > 0
                            UR = U(:,1:NR);
                            SR = S(1:NR,1:NR);
                            VR = V(:,1:NR);
                            sigma_n = sigma_n(1:NR,:);
                            R = bsxfun(@rdivide,FM*VR,(sigma_n'*Snorm));
                        end
                        
                        % Calculate CBFP weights for all requested input points (not just at the CBFP support points)
                        W = pinv(R)*FM;
                        
                        % Handle basis functions as FarField objects
                        x = FFobj.x;
                        y = FFobj.y;
                        freq = 1;
                        Prad = ones(size(freq)).*4*pi;
                        radEff = ones(size(freq));
                        coorSys = FFobj.coorSys;
                        polType = FFobj.polType;
                        gridType = FFobj.gridType;
                        freqUnit = FFobj.freqUnit;
                        freqRange = FFobj.freq;
                        
                        for ii = 1:NR
                            basis_E1 = R(1:Nang,ii);
                            basis_E2 = R(Nang+1:end,ii);
                            basis_E3 = ones([Nang,1]);
                            
                            FFbasis = FarField(x,y,basis_E1,basis_E2,basis_E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                            
                            basis(1,ii) = FFbasis;
                        end
                    else
                        %% FREQUENCY AND GEOMETRIC VARIATION
                        
                        isFreqBasis = 0;
                        
                        FFobj = reshape(FFobj,[1,nFFobjs]);
                        
                        % Extract useful parameters
                        Nf = FFobj(1,1).Nf;
                        Nang = FFobj(1,1).Nang;
                        freqHz = FFobj(1,1).freqHz;
                        
                        for ii = 1:nFFobjs
                            assert(Nf == FFobj(1,ii).Nf,'All FFobj.Nf should be identical')
                            assert(Nang == FFobj(1,ii).Nang,'All FFobj.Nang should be identical')
                        end
                        
                        %% Frequency treated as a parameter
                        
                        if Nf == 1
                            % Determine FF indices from which CBFP matrix is built                           
                            if isfield(i_FM,'x')
                                range_x = reshape(i_FM.x,1,length(i_FM.x));
                            else
                                range_x = 1:nFFobjs;
                            end
                            
                            % Get CBFP matrix
                            for ii = 1:nFFobjs
                                E1 = FFobj(1,ii).E1;
                                E2 = FFobj(1,ii).E2;
                                
                                FM_E1(:,ii) = E1;
                                FM_E2(:,ii) = E2;
                            end
                            FM = [FM_E1;FM_E2];
                            FM = FM(:,range_x);
                            
                            % Get the SVD of the FM matrix
                            [U,S,V] = svd(FM,'econ');
                        
                            Snorm = S(1,1);
                            sigma_n = diag(S)./Snorm;   % Normalization factor for sigma - this leaves an option to zero one of the field components if required and get rid of noise
                            
                            % Sort out the SVD matrix sizes due to reduced number of significant singular values
                            if nBasis < inf
                                NR = nBasis;
                            else
                                NR = length(find(sigma_n > tol));
                            end
                            
                            if NR > 0
                                UR = U(:,1:NR);
                                SR = S(1:NR,1:NR);
                                VR = V(:,1:NR);
                                sigma_n = sigma_n(1:NR,:);
                                R = bsxfun(@rdivide,FM*VR,(sigma_n'*Snorm));
                            end
                            
                            % Calculate CBFP weights for all requested input points (not just at the CBFP support points)
                            W = pinv(R)*FM;
                            
                            % Handle basis functions as FarField objects
                            x = FFobj(1,1).x;
                            y = FFobj(1,1).y;
                            freq = 1;
                            Prad = ones(size(freq)).*4*pi;
                            radEff = ones(size(freq));
                            coorSys = FFobj(1,1).coorSys;
                            polType = FFobj(1,1).polType;
                            gridType = FFobj(1,1).gridType;
                            freqUnit = FFobj(1,1).freqUnit;
                            freqRange = 1;
                            
                            for ii = 1:NR
                                basis_E1 = R(1:Nang,ii);
                                basis_E2 = R(Nang+1:end,ii);
                                basis_E3 = ones([Nang,1]);
                                
                                FFbasis = FarField(x,y,basis_E1,basis_E2,basis_E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                                
                                basis(1,ii) = FFbasis;
                            end
                            
                        else
                            %% Geometric variation per Frequency
                            
                            for ii = 1:nFFobjs
                                assert(isequal(freqHz,FFobj(1,ii).freqHz),'All FFobj.freqHz should be identical')
                            end
                            
                            % Determine FF indices from which CBFP matrix is built
                            if isfield(i_FM,'f'), range_f = reshape(i_FM.f,1,length(i_FM.f)); else, range_f = 1:Nf; end
                            if isfield(i_FM,'x'), range_x = reshape(i_FM.x,1,length(i_FM.x)); else, range_x = 1:nFFobjs; end
                            
                            c = 1;
                            for ii = range_f
                                for jj = range_x
                                    E1 = FFobj(1,jj).E1(:,ii);
                                    E2 = FFobj(1,jj).E2(:,ii);
                                    
                                    FM_E1(:,jj) = E1;
                                    FM_E2(:,jj) = E2;
                                end
                                
                                FM = [FM_E1;FM_E2];
                                
                                % Get the SVD of the FM matrix
                                [U,S,V] = svd(FM,'econ');
                                
                                Snorm = S(1,1);
                                sigma_ns = diag(S)./Snorm;   % Normalization factor for sigma - this leaves an option to zero one of the field components if required and get rid of noise
                                
                                
                                % Sort out the SVD matrix sizes due to reduced number of significant singular values
                                if nBasis < inf
                                    NR = nBasis;
                                else
                                    NR = length(find(sigma_ns > tol));
                                end
                                
                                if NR > 0
                                    URs = U(:,1:NR);
                                    SRs = S(1:NR,1:NR);
                                    VRs = V(:,1:NR);
                                    sigma_ns = sigma_ns(1:NR,:);
                                    
                                    R = bsxfun(@rdivide,FM*VRs,(sigma_ns'*Snorm));
                                    
                                    % For each frequency, compile the U, S, V matrices
                                    UR(:,:,c) = URs;
                                    SR(:,:,c) = SRs;
                                    VR(:,:,c) = VRs;
                                    sigma_n(:,:,c) = reshape(sigma_ns,[1 length(sigma_ns)]);
                                end
                                
                                % Calculate CBFP weights for all requested input points (not just at the CBFP support points)
                                Ws = pinv(R)*FM;
                                W(:,:,c) = Ws;
                                
                                % Handle basis functions as FarField objects
                                x = FFobj(1,1).x;
                                y = FFobj(1,1).y;
                                freq = FFobj(1,1).freq(1,ii);
                                Prad = ones(size(freq)).*4*pi;
                                radEff = ones(size(freq));
                                coorSys = FFobj(1,1).coorSys;
                                polType = FFobj(1,1).polType;
                                gridType = FFobj(1,1).gridType;
                                freqUnit = FFobj(1,1).freqUnit;
                                freqRange = FFobj(1,1).freq;
                                
                                for kk = 1:NR
                                    basis_E1 = R(1:Nang,kk);
                                    basis_E2 = R(Nang+1:end,kk);
                                    basis_E3 = ones([Nang,1]);
                                    
                                    FFbasis = FarField(x,y,basis_E1,basis_E2,basis_E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                                    
                                    basisFn(1,kk) = FFbasis;
                                end
                                
                                basis(:,:,c) = basisFn;
                                c = c+1;
                            end
                        end
                    end

                    %% BUILD THE CBFP OBJECT
                    obj.nBasis = NR;
                    obj.nCoeffs = NR;
                    obj.basis = basis;
                    obj.coeffs = W;
                    obj.tol = tol;
                    obj.i_FM = i_FM;
                    obj.UR = UR;
                    obj.SR = SR;
                    obj.VR = VR;
                    obj.sigma_n = sigma_n;
                    obj.isFreqBasis = isFreqBasis;
                    obj.freqRange = freqRange;
                end
            end
        end
        
        methods (Static)
            function FFobj = expansion2FarField(CBFPobj,Win,tol,i_B,nBasis)
                % function FFobj = expansion2FarField(CBFPobj,Win,tol,i_FM,nBasis)
                %
                % Inputs:   CBFPobj - CBFP object
                %           W       - Interpolated weights
                %           i_B    - Array of indices for basis functions
                %           nBasis  - Number of basis functions
                %
                % Output:   FFobj - FarField object
                
                % INPUT CHECKS
                NB = CBFPobj.nBasis;
                                
                if nargin < 2
                    Win = CBFPobj.coeffs;
                    NR = NB;
                    range_B = [1:NR];
                    
                    tol = 1e-100;
                    i_B = [];
                    nBasis = inf;
                elseif nargin < 3
                    NR = NB;
                    range_B = [1:NR];
                                        
                    tol = 1e-100;
                    i_B = [];
                    nBasis = inf;
                elseif nargin < 4
                    sigma = CBFPobj.sigma_n;
                    NR = length(find(sigma > tol));
                    range_B = [1:NR];
                                     
                    i_B = [];
                    nBasis = inf;
                elseif nargin < 5
                    assert(length(i_B) <= NB,'Number of indices specified exceeds the number of basis functions')
                    assert(max(i_B) <= NB,'Maximum index specified exceeds the number of basis functions')
                    
                    NR = length(i_B);
                    range_B = i_B;
                                    
                    nBasis = inf;
                elseif nargin < 6
                    assert(nBasis <= NB,'Specified basis functions more than the total number of basis functions')
                    NR = nBasis;
                    range_B = [1:NR];
                end
                              
                % Ignore some function inputs
                if isempty(Win), Win = CBFPobj.coeffs; freq = CBFPobj.freqRange; end
                if isempty(tol), tol = 1e-100; end
                if isempty(i_B), i_B = []; end
                if isempty(nBasis), nBasis = inf; end
                
                freqBasis = CBFPobj.isFreqBasis;
                                
                B_E1 = [];
                B_E2 = [];
                
                % Expansion across Frequency
                if freqBasis
                    nCoeffs = length(Win(1,:));
                    
                    if nargin < 2,  freq = CBFPobj.freqRange; else, if ~isempty(Win), freq = ones(1,nCoeffs); end, end
                    
                    for ii = range_B
                        % Get the Basis Functions
                        basisFn = CBFPobj.basis(1,ii);
                        
                        B_E1 = [B_E1 basisFn.E1];
                        B_E2 = [B_E2 basisFn.E2];
                    end
                    B_E = [B_E1;B_E2];
                    
                    FF = B_E * Win(range_B,:);
                    
                    % Build FarField Object
                    Nang = basisFn.Nang;
                    
                    x = basisFn.x;
                    y = basisFn.y;
                    
                    E1 = FF(1:Nang,:);
                    E2 = FF(Nang+1:end,:);
                    E3 = zeros(size(E1));
                    
                    Prad = ones([1,nCoeffs]).*4*pi;
                    radEff = ones([1,nCoeffs]);
                    coorSys = basisFn.coorSys;
                    polType = basisFn.polType;
                    gridType = basisFn.gridType;
                    freqUnit = basisFn.freqUnit;
                    
                    
                    FFobj = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                    
%                     FFobj = setEnames(FFobj);
%                     FFobj = setXYnames(FFobj);
%                     FFobj = setPhTh(FFobj);
%                     FFobj = setFreq(FFobj);
%                     FFobj = setBase(FFobj);
                else
                    dim = length(size(CBFPobj.basis));
                    
                    switch dim
                        case 2
                            % Frequency treated as a parameter
                            nCoeffs = length(Win(1,:));
                            
                            for ii = range_B
                                % Get the Basis Functions
                                basisFn = CBFPobj.basis(1,ii);
                                
                                B_E1 = [B_E1 basisFn.E1];
                                B_E2 = [B_E2 basisFn.E2];
                            end
                            B_E = [B_E1;B_E2];
                            FF = B_E * Win(range_B,:);
                            
                            % Build FarField Objects
                            Nang = basisFn.Nang;
                            
                            x = basisFn.x;
                            y = basisFn.y;
                                                       
                            freq = 1;
                            Prad = 4*pi;
                            radEff = 1;
                            coorSys = basisFn.coorSys;
                            polType = basisFn.polType;
                            gridType = basisFn.gridType;
                            freqUnit = basisFn.freqUnit;
                            
                            for ii = 1:nCoeffs
                                E1 = FF(1:Nang,ii);
                                E2 = FF(Nang+1:end,ii);
                                E3 = zeros(size(E1));
                                FFobj(1,ii) = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                            end
                            
                        case 3
                            % Geometric variation
                            nCoeffs = length(Win(1,:,1));
                            Nf = length(CBFPobj.basis(1,1,:));
                            
                            % Build FarField Objects
                            Nang = CBFPobj.basis(1,1,1).Nang;
                            x = CBFPobj.basis(1,1,1).x;
                            y = CBFPobj.basis(1,1,1).y;
                            
                            Prad = 4*pi;
                            radEff = 1;
                            coorSys = CBFPobj.basis(1,1,1).coorSys;
                            polType = CBFPobj.basis(1,1,1).polType;
                            gridType = CBFPobj.basis(1,1,1).gridType;
                            freqUnit = CBFPobj.basis(1,1,1).freqUnit;
                            
                            for ii = 1:Nf
                                B_E1 = [];
                                B_E2 = [];
                                B_E = [];
                                for jj = range_B
                                    % Get the Basis Functions
                                    basisFn = CBFPobj.basis(1,jj,ii);
                                    
                                    B_E1 = [B_E1 basisFn.E1];
                                    B_E2 = [B_E2 basisFn.E2];
                                end
                                B_E = [B_E1;B_E2];
                                FF = B_E * Win(range_B,:,ii);
                                
                                freq = CBFPobj.freqRange(1,ii);
                                
                                for kk = 1:nCoeffs
                                    E1 = FF(1:Nang,kk);
                                    E2 = FF(Nang+1:end,kk);
                                    E3 = zeros(size(E1));
                                    FFobjs(1,kk) = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                                end
                                
                                FFobj(1,:,ii) = FFobjs;
                            end
                    end
                end
            end
        end
end