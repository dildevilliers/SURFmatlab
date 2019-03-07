classdef CBFP < FarFieldExpansion 
    properties
        tol
        i_FM 
        UR
        S
        VR
        sigma_n
    end
    
        methods
            function obj = CBFP(FFobj,tol,i_FM,nBasis)
                % function obj = CBFP(FFobj,tol,i_FM,nBasis)
                % CBFP object constructor method
                %
                % Inputs:   FFobj   - Farfield object
                %           tol     - Tolerance 
                %           i_FM    - SVD: reduced left singular vectors
                %           nBasis  - Number of basis functions
                %
                % Output:   CBFPobj - CBFP object
                
                if nargin == 0
                    % Expand gausian beam pattern
                else
                    % CHECKS
                    
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
                            VR = V(:,1:NR);
                            R = bsxfun(@rdivide,FM*VR,(sigma_n(1:NR)'*Snorm));
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
                        
                        for ii = 1:NR
                            basis_E1 = R(1:Nang,ii);
                            basis_E2 = R(Nang+1:end,ii);
                            basis_E3 = ones([Nang,1]);
                            
                            FFbasis = FarField(x,y,basis_E1,basis_E2,basis_E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                            
                            basis(1,ii) = FFbasis;
                        end
                    else
                        %% FREQUENCY AND GEOMETRIC VARIATION
                        
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
                                VR = V(:,1:NR);
                                R = bsxfun(@rdivide,FM*VR,(sigma_n(1:NR)'*Snorm));
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
                                sigma_n(:,:,c) = reshape(sigma_ns,[1 length(sigma_ns)]);
                                
                                % Sort out the SVD matrix sizes due to reduced number of significant singular values
                                if nBasis < inf
                                    NR = nBasis;
                                else
                                    NR = length(find(sigma_ns > tol));
                                end
                                
                                if NR > 0
                                    URs = U(:,1:NR);
                                    VRs = V(:,1:NR);
                                    R = bsxfun(@rdivide,FM*VRs,(sigma_ns(1:NR)'*Snorm));
                                    
                                    % For each frequency, compile the U, S, V matrices
                                    UR(:,:,c) = URs;
                                    S(:,:,c) = S;
                                    VR(:,:,c) = VRs;
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
                    obj.S = S;
                    obj.VR = VR;
                    obj.sigma_n = sigma_n;
                end
            end
            
%             function [FFobj] = Expansion2FarField(CBFPobj,W,i_FM,nBasis) 
            
        end
end