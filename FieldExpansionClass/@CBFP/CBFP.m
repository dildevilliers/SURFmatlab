classdef CBFP < FarFieldExpansion
    
    % Class: CBFP (Characteristic Basis Function Patterns)
    % Superclass: FarFieldExpansion (abstract)
    % Subclass(es): none
    % Methods allow transformation of a target far-field(s) to a set of basis functions
    % and corresponding coefficients 
    % Farfields can be reconstructed from the generated basis functions and coefficients 
    
    properties
        tol
        iFM 
    end
    
    properties (SetAccess = private)
        nBasis
        basis
        nCoeffs
        coeffs
        UR
        SR
        VR
        sigma_n
    end
    
    properties (SetAccess = immutable, Hidden = true)
        flagMode
        freqRange
    end
    
    methods
        function CBFPobj = CBFP(FFobj,tol,iFM,iBasis)
            % function obj = CBFP(FFobj,tol,iFM,iBasis)
            %
            % CBFP constructor method
            % Pattern variation: Frequency - Parse a farfield object defined at multiple frequencies
            %                    Geometric - Parse an array of [1xN] farfield objects
            %                                Single Freq objects: Assumption is frequency is treated as a normal paramter
            %                                Multiple Freq objects: Frequencies should be the same at all objects. 
            %                                                       Expansion is done per frequency for geometric variation
            % Inputs:   
            %-- FFobj:  Farfield object containing a set of farfields
            %           [1xN] FFobjs or a single FFobj
            %-- tol:    Max singular value tolerance of the basis functions to keep
            %           Default 1e-100
            %-- iFM:    Struct specifying which farfield patterns to use in generating basis functions
            %           iFM.f => [1xN] indices selecting farfields across the frequency axis 
            %           iFM.x => [1xN] indices selecting farfields across the design space
            %-- iBasis: Indices to select which basis functions to return in the CBFP object
            %           [1xN] row vector 
            %
            % Output:   CBFPobj - CBFP object
            
            if nargin == 0
                % Expand gausian beam pattern
            else
                % INPUT CHECKS
                if nargin < 2
                    tol = 1e-100;
                    iFM = [];
                    iBasis = [];
                elseif nargin < 3
                    iFM = [];
                    iBasis = [];
                elseif nargin < 4
                    iBasis = [];
                end
                
                % Ignore some function inputs
                if isempty(tol), tol = 1e-100; end
                if isempty(iFM), iFM = []; end
                if isempty(iBasis), iBasis = []; end
                
                nFFobjs = length(FFobj);
                if nFFobjs == 1
                    if isfield(iFM,'x'), warning('No geometric variation. iFM.x ignored'), end
                    if isfield(iFM,'f'),assert(length(iFM.f) <= FFobj.Nf,'Specified number of indices,iFM.f larger than number of frequency samples,FFobj.Nf'), end
                    if isfield(iFM,'f'),assert(max(iFM.f) <= FFobj.Nf,'Largest specified index in iFM.f larger than maximum possible index FFobj.Nf'), end
                    if ~isempty(iBasis), assert(max(iBasis)<=FFobj.Nf,'Maximum index specified exceeds the largest number of basis functions that can be generated'); end
                else
                    
                    Nf = FFobj(1,1).Nf;
                    if Nf == 1
                        if isfield(iFM,'f'), warning('Frequency treated as a normal parameter. iFM.f ignored'), end
                        if isfield(iFM,'x'),assert(length(iFM.x) <= nFFobjs,'Specified number of indices,iFM.x larger than number of geometric variations, num of FFobjs'), end
                        if isfield(iFM,'x'),assert(max(iFM.x) <= nFFobjs,'Largest specified index in iFM.x larger than maximum possible index, num of FFobjs'), end
                        if ~isempty(iBasis), assert(max(iBasis)<=nFFobjs,'Maximum index specified exceeds the largest number of basis functions that can be generated'); end
                    else
                        if isfield(iFM,'f'),assert(length(iFM.f) <= Nf,'Specified number of indices,iFM.f larger than number of frequency samples,FFobj.Nf'), end
                        if isfield(iFM,'f'),assert(max(iFM.f) <= Nf,'Largest specified index in iFM.f larger than maximum possible index FFobj.Nf'), end
                        
                        if isfield(iFM,'x'),assert(length(iFM.x) <= nFFobjs,'Specified number of indices,iFM.x larger than number of geometric variations, num of FFobjs'), end
                        if isfield(iFM,'x'),assert(max(iFM.x) <= nFFobjs,'Largest specified index in iFM.x larger than maximum possible index, num of FFobjs'), end
                        if ~isempty(iBasis), assert(max(iBasis)<=nFFobjs,'Maximum index specified exceeds the largest number of basis functions that can be generated'); end
                    end
                end
                
                %% FREQUENCY ONLY VARIATION
        
                if nFFobjs == 1
                    mode=1; % flag - freq variation
                    
                    % Extract useful parameters
                    Nf = FFobj.Nf;
                    E1 = FFobj.E1;
                    E2 = FFobj.E2;
                    Nang = FFobj.Nang;
                    
                    % Determine FF indices from which CBFP matrix is built
                    if isfield(iFM,'f')
                        range_f = reshape(iFM.f,1,length(iFM.f));
                    else
                        range_f = 1:Nf;
                    end                 
                    % Get the FM matrix
                    FM_E1 = [];
                    FM_E2 = [];
                    FM_E1 = [FM_E1 E1(:,range_f)];
                    FM_E2 = [FM_E2 E2(:,range_f)];
                    
                    FM = [FM_E1;FM_E2];
                    
                    % Get the SVD of the FM matrix
                    [U,S,V] = svd(FM,'econ');
                    
                    Snorm = S(1,1);
                    sigma_n_mat = diag(S)./Snorm;   % Normalization factor for sigma - this leaves an option to zero one of the field components if required and get rid of noise
                    
                    % Sort out the SVD matrix sizes due to reduced number of significant singular values
                    if ~isempty(iBasis)
                        NR = length(iBasis);
                        range_B = iBasis;
                    else
                        NR = length(find(sigma_n_mat > tol));
                        range_B = (1:NR);
                    end
                    
                    if NR > 0
                        Umat = U(:,range_B);
                        Smat = S(range_B,range_B);
                        Vmat = V(:,range_B);
                        sigma_n_mat = sigma_n_mat(range_B,:);
                    end
                    
                    % Calculate CBFP weights for all requested input points (not just at the CBFP support points)
                    W{1} = pinv(Umat)*FM;
                    UR{1} = Umat;
                    SR{1} = Smat;
                    VR{1} = Vmat;
                    sigma_n{1} = sigma_n_mat;
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
                    freqRange = FFobj.freqHz;
                    
                    for ii = 1:NR
                        basis_E1 = Umat(1:Nang,ii);
                        basis_E2 = Umat(Nang+1:end,ii);
                        basis_E3 = ones([Nang,1]);
                        
                        FFbasis = FarField(x,y,basis_E1,basis_E2,basis_E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                        
                        basis{1,ii} = FFbasis;
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
                        mode=2; % flag - freq treated as a parameter
                        
                        % Determine FF indices from which CBFP matrix is built
                        if isfield(iFM,'x')
                            range_x = reshape(iFM.x,1,length(iFM.x));
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
                        sigma_n_mat = diag(S)./Snorm;   % Normalization factor for sigma - this leaves an option to zero one of the field components if required and get rid of noise
                        
                        % Sort out the SVD matrix sizes due to reduced number of significant singular values
                        if ~isempty(iBasis)
                            NR = length(iBasis);
                            range_B = iBasis;
                        else
                            NR = length(find(sigma_n_mat > tol));
                            range_B = (1:NR);
                        end
                        if NR > 0
                            Umat = U(:,range_B);
                            Smat = S(range_B,range_B);
                            Vmat = V(:,range_B);
                            sigma_n_mat = sigma_n_mat(range_B,:);
                        end
                        
                        % Calculate CBFP weights for all requested input points (not just at the CBFP support points)
                        Wmat = pinv(Umat)*FM;
                        W = mat2cell(Wmat,length(Wmat(:,1)),ones(1,length(Wmat(1,:))));
                        UR{1} = Umat;
                        SR{1} = Smat;
                        VR{1} = Vmat;
                        sigma_n{1} = sigma_n_mat;
                        
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
                            basis_E1 = Umat(1:Nang,ii);
                            basis_E2 = Umat(Nang+1:end,ii);
                            basis_E3 = ones([Nang,1]);
                            
                            FFbasis = FarField(x,y,basis_E1,basis_E2,basis_E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                            
                            basis{1,ii} = FFbasis;
                        end
                        
                    else
                        %% Geometric variation per Frequency
                        mode = 3; % flag - per freq geometric variation
                        
                        for ii = 1:nFFobjs
                            assert(isequal(freqHz,FFobj(1,ii).freqHz),'All FFobj.freq should be identical')
                        end
                        
                        % Determine FF indices from which CBFP matrix is built
                        if isfield(iFM,'f'), range_f = reshape(iFM.f,1,length(iFM.f)); else, range_f = 1:Nf; end
                        if isfield(iFM,'x'), range_x = reshape(iFM.x,1,length(iFM.x)); else, range_x = 1:nFFobjs; end
                        
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
                            if ~isempty(iBasis)
                                NR = length(iBasis);
                                range_B = iBasis;
                            else
                                NR = length(find(sigma_ns > tol));
                                range_B = (1:NR);
                            end
                            if NR > 0
                                URs = U(:,range_B);
                                SRs = S(range_B,range_B);
                                VRs = V(:,range_B);
                                sigma_ns = sigma_ns(range_B,:);
                                
                                % For each frequency, compile the U, S, V matrices
                                URmat(:,:,c) = URs;
                                SRmat(:,:,c) = SRs;
                                VRmat(:,:,c) = VRs;
                                sigma_n_mat(:,:,c) = reshape(sigma_ns,[1 length(sigma_ns)]);
                            end
                            UR{1} = URmat;
                            SR{1} = SRmat;
                            VR{1} = VRmat;
                            sigma_n{1} = sigma_n_mat;
                            
                            % Calculate CBFP weights for all requested input points (not just at the CBFP support points)
                            Ws = pinv(URs)*FM;
                            b=1;
                            for jj = range_x
                                W{1,jj}(:,c) = Ws(:,b);
                                b=b+1;
                            end
                                                        
                            % Handle basis functions as FarField objects
                            x = FFobj(1,1).x;
                            y = FFobj(1,1).y;
                            freq = FFobj(1,1).freqHz(1,ii);
                            Prad = ones(size(freq)).*4*pi;
                            radEff = ones(size(freq));
                            coorSys = FFobj(1,1).coorSys;
                            polType = FFobj(1,1).polType;
                            gridType = FFobj(1,1).gridType;
                            freqUnit = FFobj(1,1).freqUnit;
                            freqRange = FFobj(1,1).freqHz;
                            
                            for kk = 1:NR
                                basis_E1 = URs(1:Nang,kk);
                                basis_E2 = URs(Nang+1:end,kk);
                                basis_E3 = ones([Nang,1]);
                                
                                FFbasis = FarField(x,y,basis_E1,basis_E2,basis_E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                                
                                basis{c,kk} = FFbasis;
                            end
                            c = c+1;
                        end
                    end
                end
                
                % BUILD THE CBFP OBJECT
                CBFPobj.nBasis = NR;
                CBFPobj.nCoeffs = NR;
                CBFPobj.basis = basis;
                CBFPobj.coeffs = W;
                CBFPobj.tol = tol;
                CBFPobj.iFM = iFM;
                CBFPobj.UR = UR;
                CBFPobj.SR = SR;
                CBFPobj.VR = VR;
                CBFPobj.sigma_n = sigma_n;
                CBFPobj.flagMode = mode;
                CBFPobj.freqRange = freqRange;
            end
        end
    end
    
    methods (Static)
        function Wout = farField2Coeffs(FFobj,CBFPobj,tol,iBasis)
            % Wout = farField2Coeffs(FFobj,CBFPobj,tol,iBasis)
            %
            % Target field coefficients from a pre-exisiting set of basis functions
            %
            % Inputs:   
            %-- FFobj:      Farfield object containing a set of target farfields
            %               [1xN] FFobjs or a single FFobj
            %-- CBFPobj:    Previously constructed set of basis functions
            %-- tol:        Max singular value tolerance of the basis functions to keep
            %               Default 1e-100
            %-- iBasis:     Indices to select which basis functions to rebuild the target field with
            %               [1xN] row vector 
            %
            % Output:   	Wout - A set of coefficients for each target farfield
            %               {1xN} cell array or {NfxN} cell array - Only when expanding farfields for geometric variation per frequency

            % INPUT CHECKS
            NB = CBFPobj.nBasis;

            if nargin < 2
                error('Target field (FFobj) and basis functions (CBFPobj) required');           
            elseif nargin < 3
                range_B = (1:NB);
                tol = 1e-100;
                iBasis = [];
            elseif nargin < 4
                sigma = CBFPobj.sigma_n{:,1};
                range_B = (1:length(find(sigma > tol)));
                iBasis = [];
            elseif nargin < 5
                assert(length(iBasis) <= NB,'Number of indices specified exceeds the number of basis functions')
                range_B = iBasis;
            end
            
            % Ignore some function inputs            
            if isempty(tol), tol = 1e-100; end
            if isempty(iBasis), iBasis = []; range_B = (1:NB); else, assert(max(iBasis) <= NB,'Maximum index specified exceeds the number of basis functions'), end
            
            nFFobjs = length(FFobj);
            mode = CBFPobj.flagMode;
            freqRnge = CBFPobj.freqRange;
            switch mode
                case 1
                    if nFFobjs == 1
                        freqs = FFobj.freqHz;
                        Nf = FFobj.Nf;
                        
                        if (max(freqRnge) < max(freqs)), warning('Specified maximum frequency of the target field is outside the frequency range that the basis functions were constructed from'), end
                        
                        B_E1 = [];
                        B_E2 = [];
                        for ii = range_B
                            % Get the Basis Functions
                            basisFn = CBFPobj.basis{1,ii};
                            
                            B_E1 = [B_E1 basisFn.E1];
                            B_E2 = [B_E2 basisFn.E2];
                        end
                        B_E = [B_E1;B_E2];
                        
                        for ii = 1:Nf
                            E1(:,ii) = FFobj.E1(:,ii);
                            E2(:,ii) = FFobj.E2(:,ii);
                        end
                        FF = [E1;E2];  
                        Wmat = pinv(B_E)*FF;
                        Wout{1} = Wmat;
                    else
                        B_E1 = [];
                        B_E2 = [];
                        for ii = range_B
                            % Get the Basis Functions
                            basisFn = CBFPobj.basis{1,ii};
                            
                            B_E1 = [B_E1 basisFn.E1];
                            B_E2 = [B_E2 basisFn.E2];
                        end
                        B_E = [B_E1;B_E2];
                        
                        for ii = 1:nFFobjs
                            FF = FFobj(1,ii);
                            assert(FF.Nf <= 1,'Multiple frequencies defined. Specify one frequency at a time')
                            
                            if (max(freqRnge) < max(FF.freq)), warning('Specified maximum frequency of the target field is outside the frequency range that the basis functions were constructed from'), end
                            
                            E1 = FF.E1(:,1);
                            E2 = FF.E2(:,1);
                            FF_E = [E1;E2];
                            
                            Wout{1,ii} = pinv(B_E)*FF_E;
                        end
                    end
                case 2
                    if nFFobjs == 1
                        assert(FFobj.Nf <= 1,'Multiple frequencies defined. Specify one frequency at a time')
                        
                        B_E1 = [];
                        B_E2 = [];
                        for ii = range_B
                            % Get the Basis Functions
                            basisFn = CBFPobj.basis{1,ii};
                            
                            B_E1 = [B_E1 basisFn.E1];
                            B_E2 = [B_E2 basisFn.E2];
                        end
                        B_E = [B_E1;B_E2];
                        
                        E1 = FFobj.E1;
                        E2 = FFobj.E2;
                        FF = [E1;E2];
                        
                        Wmat = pinv(B_E)*FF;
                        Wout{1} = Wmat;
                    else
                        B_E1 = [];
                        B_E2 = [];
                        for ii = range_B
                            % Get the Basis Functions
                            basisFn = CBFPobj.basis{1,ii};
                            
                            B_E1 = [B_E1 basisFn.E1];
                            B_E2 = [B_E2 basisFn.E2];
                        end
                        B_E = [B_E1;B_E2];
                        
                        for ii = 1:nFFobjs
                            FF = FFobj(1,ii);
                            assert(FF.Nf <= 1,'Multiple frequencies defined. Specify one frequency at a time in the FarField object')
                            
                            E1 = FF.E1(:,1);
                            E2 = FF.E2(:,1);
                            FF_E = [E1;E2];
                            Wmat = pinv(B_E)*FF_E;
                            
                            Wout{1,ii} = Wmat;
                        end
                    end
                case 3
                    if nFFobjs == 1
                        Nf = FFobj.Nf;
                        
                        if Nf == 1
                            k = ismember(freqRnge,FFobj.freqHz);
                            idx = find(k,1);
                            
                            assert(~isempty(idx),'Frequency of the Farfield object does not correspond to any of the frequencies used to build the basis functions')
                            
                            basisFns = CBFPobj.basis(idx,:);
                            
                            B_E1 = [];
                            B_E2 = [];
                            for ii = range_B
                                % Get the Basis Functions
                                basisFn = basisFns{1,ii};
                                
                                B_E1 = [B_E1 basisFn.E1];
                                B_E2 = [B_E2 basisFn.E2];
                            end
                            B_E = [B_E1;B_E2];
                            
                            E1 = FFobj.E1;
                            E2 = FFobj.E2;
                            FF = [E1;E2];
                            
                            Wmat = pinv(B_E)*FF;
                            Wout{1} = Wmat;
                        else
                            for ii = 1:Nf
                                k = ismember(freqRnge,FFobj.freqHz(1,ii));
                                idx = find(k,1);
                                
                                assert(~isempty(idx),'Frequency of the Farfield object does not correspond to any of the frequencies used to build the basis functions')
                                
                                basisFns = CBFPobj.basis(idx,:);
                                
                                B_E1 = [];
                                B_E2 = [];
                                for jj = range_B
                                    % Get the Basis Functions
                                    basisFn = basisFns{1,jj};
                                    
                                    B_E1 = [B_E1 basisFn.E1];
                                    B_E2 = [B_E2 basisFn.E2];
                                end
                                B_E = [B_E1;B_E2];
                                
                                E1 = FFobj.E1(:,ii);
                                E2 = FFobj.E2(:,ii);
                                FF = [E1;E2];
                                
                                Wmat(:,ii) = pinv(B_E)*FF;
                            end
                            Wout{1} = Wmat;
                        end
                    else
                        for ii = 1:nFFobjs
                            FFobj_1 = FFobj(1,ii);
                            Nf = FFobj_1.Nf;
                            Wmat = [];
                            
                            if Nf == 1
                                k = ismember(freqRnge,FFobj_1.freq);
                                idx = find(k,1);
                                
                                assert(~isempty(idx),'Frequency of the Farfield object does not correspond to any of the frequencies used to build the basis functions')
                                
                                basisFns = CBFPobj.basis(idx,:);
                                
                                B_E1 = [];
                                B_E2 = [];
                                for jj = range_B
                                    % Get the Basis Functions
                                    basisFn = basisFns{1,jj};
                                    
                                    B_E1 = [B_E1 basisFn.E1];
                                    B_E2 = [B_E2 basisFn.E2];
                                end
                                B_E = [B_E1;B_E2];
                                
                                E1 = FFobj.E1;
                                E2 = FFobj.E2;
                                FF = [E1;E2];
                                
                                Wmat = pinv(B_E)*FF;
                            else
                                for jj = 1:Nf
                                    k = ismember(freqRnge,FFobj_1.freqHz(1,jj));
                                    idx = find(k,1);
                                    
                                    assert(~isempty(idx),'Frequency of the Farfield object does not correspond to any of the frequencies used to build the basis functions')
                                    
                                    basisFns = CBFPobj.basis(idx,:);
                                    
                                    B_E1 = [];
                                    B_E2 = [];
                                    for kk = range_B
                                        % Get the Basis Functions
                                        basisFn = basisFns{1,kk};
                                        
                                        B_E1 = [B_E1 basisFn.E1];
                                        B_E2 = [B_E2 basisFn.E2];
                                    end
                                    B_E = [B_E1;B_E2];
                                    
                                    E1 = FFobj_1.E1(:,jj);
                                    E2 = FFobj_1.E2(:,jj);
                                    FF = [E1;E2];
                                    
                                    Wmat(:,jj) = pinv(B_E)*FF;
                                end
                            end
                            Wout{1,ii} = Wmat;
                        end
                    end
            end
        end
        function FFobj = coeffs2FarField(CBFPobj,W,tol,iBasis,freqIndex)
            % function FFobj = coeffs2FarField(CBFPobj,W,tol,iBasis,freqIndex)
            %
            % Reconstruct FarFields from previously generated basis functions
            % Inputs:   
            %-- CBFPobj:    Previously constructed set of basis functions
            %-- W:          Basis function coefficients. These may be from an interpolation function
            %               {1xN} cell array or {NfxN} cell array - Only when expanding farfields for geometric variation per frequency
            %-- tol:        Max singular value tolerance of the basis functions to keep
            %               Default 1e-100
            %-- iBasis:     Indices to select which basis functions to rebuild the FarField object with
            %               [1xN] row vector 
            %-- freqIndex:  Optional input required and used only when CBFP expansion was done for geometric variation per frequency
            %               Specifies which frequency slice to reconstruct a farfield pattern at 
            %               To be specified when coefficents parsed to the method do not cover the total number of frequency points across the frequency axis
            %               [1xN] integer array  
            %
            %-- FFobj:  Farfield object
            
            % INPUT CHECKS
            NB = CBFPobj.nBasis;
            if nargin < 2
                W = CBFPobj.coeffs;
                range_B = (1:NB);
                tol = 1e-100;
                iBasis = [];
                freqIndex = [];
            elseif nargin < 3
                range_B = (1:NB);
                tol = 1e-100;
                iBasis = [];
                freqIndex = [];
            elseif nargin < 4
                sigma = CBFPobj.sigma_n{:,1};
                range_B = (1:length(find(sigma > tol)));
                iBasis = [];
                freqIndex = [];
            elseif nargin < 5
                assert(length(iBasis) <= NB,'Number of indices specified exceeds the number of basis functions')
                range_B = iBasis;
                freqIndex = [];
            end
            
            % Ignore some function inputs
            if isempty(W), W = CBFPobj.coeffs; freq = CBFPobj.freqRange; end
            if isempty(tol), tol = 1e-100; end
            if isempty(iBasis), iBasis = []; range_B = (1:NB); else, assert(max(iBasis) <= NB,'Maximum index specified exceeds the number of basis functions'), end
            if isempty(freqIndex), freqIndex = []; range_B = (1:NB); end
            
            mode = CBFPobj.flagMode;
            switch mode
                case 1
                    for ii = 1:length(W(1,:))
                        Wmat = W{1,ii};
                        nCoeffs = length(Wmat(1,:));
                        if nargin < 2,  freq = CBFPobj.freqRange; else, if ~isempty(W), freq = ones(1,nCoeffs); end, end
                        
                        assert(NB == length(Wmat(:,1)),'Column entries of W not equal to nBasis. Coefficients for all basis functions should be provided')
                        
                        B_E1 = [];
                        B_E2 = [];
                        for jj = range_B
                            % Get the Basis Functions
                            basisFn = CBFPobj.basis{1,jj};
                            
                            B_E1 = [B_E1 basisFn.E1];
                            B_E2 = [B_E2 basisFn.E2];
                        end
                        B_E = [B_E1;B_E2];
                        
                        FF = B_E * Wmat(range_B,:);
                        
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

                        FFobj(1,ii) = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                    end
                case 2
                    % Frequency treated as a parameter
                    for ii = 1:length(W(1,:))
                        Wmat = W{1,ii};
                        nCoeffs = length(Wmat(1,:));
                        assert(NB == length(Wmat(:,1)),'Column entries of W not equal to nBasis. Coefficients for all basis functions should be provided')
                        
                        B_E1 = [];
                        B_E2 = [];
                        for jj = range_B
                            % Get the Basis Functions
                            basisFn = CBFPobj.basis{1,jj};
                            
                            B_E1 = [B_E1 basisFn.E1];
                            B_E2 = [B_E2 basisFn.E2];
                        end
                        B_E = [B_E1;B_E2];
                        FF = B_E * Wmat(range_B,:);
                        
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

                        E1 = FF(1:Nang,:);
                        E2 = FF(Nang+1:end,:);
                        E3 = zeros(size(E1));
                        FFobj(1,ii) = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                    end
                case 3
                    % Geometric variation per freq
                    for ii = 1:length(W(1,:))
                        Wmat = W{1,ii};
                        nCoeffs = length(Wmat(1,:));
                        assert(NB == length(Wmat(:,1)),'Column entries of W not equal to nBasis. Coefficients for all basis functions should be provided')
                        Nf = length(CBFPobj.basis(:,1));

                        if(nCoeffs ~= Nf) && (isempty(freqIndex)), error('W column entries not equal to total number of frequencies. Specify at which frequencies to rebuild through freqIndex'), end
                        
                        if ~isempty(freqIndex)
                            assert(length(freqIndex) == nCoeffs,'Number of indices on freqIndex should be equal to the total number of columns in W')
                            range_f = freqIndex;
                        else
                            range_f = (1:Nf);
                        end
                        
                        B_E1 = [];
                        B_E2 = [];
                        B_E = [];
                        for kk = range_B
                            % Get the Basis Functions
                            basisFn = CBFPobj.basis{range_f,kk};
                            
                            B_E1 = [B_E1 basisFn.E1];
                            B_E2 = [B_E2 basisFn.E2];
                        end
                        B_E = [B_E1;B_E2];
                        FF = B_E * Wmat(range_B,:);
                        
                        % Build FarField Objects
                        Nang = CBFPobj.basis{1,1}.Nang;
                        x = CBFPobj.basis{1,1}.x;
                        y = CBFPobj.basis{1,1}.y;
                        Prad = 4*pi*ones(1,length(range_f));
                        radEff = 1*ones(1,length(range_f));
                        freq = CBFPobj.freqRange(1,range_f);
                        coorSys = CBFPobj.basis{1,1}.coorSys;
                        polType = CBFPobj.basis{1,1}.polType;
                        gridType = CBFPobj.basis{1,1}.gridType;
                        freqUnit = CBFPobj.basis{1,1}.freqUnit;
                        
                        E1(:,:) = FF(1:Nang,:);
                        E2(:,:) = FF(Nang+1:end,:);
                        E3 = zeros(size(E1));
                        FFobj(1,ii) = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
                    end
            end
        end
    end
end