classdef FarField
    % Based largely on information in CST help file: Farfield Calculation
    % Overview, and the AMTA paper in this folder 'COORDINATE SYSTEM
    % PLOTTING FOR ANTENNA MEASUREMENTS', GF Masters and SF Gregson
    % Currently assumes direction of propagation is the global z-axis for
    % all directional polarization types.
    
    properties
        r(1,1) double {mustBeReal, mustBeFinite} = 1    % Radius very E-field is evaluated in (m)
        Prad(1,:) double {mustBeReal, mustBeFinite} = 4*pi
        radEff(1,:) double {mustBeReal, mustBeFinite} = 1
        slant(1,1) double {mustBeReal, mustBeFinite} = pi/4   % slant angle in radians - measured between Exp and E1
        freqUnit(1,:) char {mustBeMember(freqUnit,{'Hz','kHz','MHz','GHz','THz'})} = 'Hz'
    end
    
    properties (SetAccess = private)
        x(:,1) double {mustBeReal, mustBeFinite}
        y(:,1) double {mustBeReal, mustBeFinite}
        E1(:,:) double %{mustBeFinite}
        E2(:,:) double %{mustBeFinite}
        E3(:,:) double %{mustBeFinite}
        freq(1,:) double {mustBeReal, mustBeFinite} = 1
        coorType(1,:) char {mustBeMember(coorType,{'spherical','Ludwig1','Ludwig2AE','Ludwig2EA','Ludwig3'})} = 'spherical'
        polType(1,:) char {mustBeMember(polType,{'linear','circular','slant'})} = 'linear'
        gridType(1,:) char {mustBeMember(gridType,{'PhTh','DirCos','AzEl','ElAz','TrueView','ArcSin'})} = 'PhTh'
        symmetryXZ = 'none'   % Type of symmetry about the xz-plane (could be electric|none|magnetic)
        symmetryYZ = 'none'   % Type of symmetry about the yz-plane (could be electric|none|magnetic)
        symmetryXY = 'none'   % Type of symmetry about the xy-plane (could be electric|none|magnetic)
        symmetryBOR1 = false % Is the pattern a BOR1 type pattern
    end
    
    properties (Dependent = true)
        ph
        th
        xname
        yname
        E1name  % ['Eth', 'Ex', 'Eaz', 'Eal', 'Eh', 'Elh', 'Exp'] - depends on coorType and polType
        E2name  % ['Eph', 'Ey', 'Eel', 'Eep', 'Ev', 'Erh', 'Eco'] - depends on coorType and polType
        Nf      % number of frequencies
        Nx     % number of unique x angles
        Ny     % number of unique y angles
        Nang    % number of angle combinations [Nx x Ny]
        freqHz
        Directivity_dBi % directivity in dBi [1 x Nf]
        Gain_dB         % Gain in dB [1 x Nf]
        radEff_dB
        xRangeType     % 'sym' or 'pos'
        yRangeType     % '180' or '360'
    end
    
    properties (SetAccess = private, Hidden = true)
        % Keep the input data here to not lose some info when going through
        % a DirCos projection and back...
        xBase
        yBase
        phBase
        thBase
        NxBase
        NyBase
        gridTypeBase
        E1Base
        E2Base
        E3Base
        coorTypeBase
        polTypeBase
        symXZ = 0
        symYZ = 0
        symXY = 0
    end
    
    properties (Constant = true, Hidden = true)
        c0 = physconst('Lightspeed');
        eps0 = 8.854187817000001e-12;
        mu0 = 1.256637061435917e-06;
        eta0 = 3.767303134749689e+02;
        nSigDig = 8;
    end
    
    methods
        % Make a basic constructor method
        function obj = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorType,polType,gridType,freqUnit,slant)
            
            % function obj = FarField(th,ph,E1,E2,E3,freq,Prad,radEff)
            % Constructor method for the FarField object
            % Required inputs:
            % x: column vector [Nang x 1] of ph angles in rad
            % y: column vector [Nang x 1] of th angles in rad
            % E1: First E-field pattern component (exp(jkr)/r suppressed) of size [Nang x Nf]
            % E2: Second E-field pattern component (exp(jkr)/r suppressed) of size [Nang x Nf]
            %
            % Optional inputs - can be left empty or omitted
            % 'E3': Third E-field pattern component (exp(jkr)/r suppressed) of size [Nang x Nf]
            % 'freq': vector [1 x Nf] of frequencies where the fields are defined in Hz
            % 'Prad': vector [1 x Nf] of radiated powers in W
            % 'radEff': vector [1 x Nf] of radiation efficiencies
            
            if nargin == 0 %No-inputs case - generate a Gaussian beam pattern at a single frequency (1 GHz) over the full sphere, 5 degree angular resolution
                Nth_cut = 37;
                Nph_cut = 73;
                th = linspace(0,pi,Nth_cut);
                ph = linspace(0,2*pi,Nph_cut);
                th0 = deg2rad(50);
                taper_dB = -10;
                freq = 1e9;
                [PH,TH] = meshgrid(ph,th);
                P = powerPattern(PH(:),TH(:),'gauss',th0,taper_dB,freq);
                obj = FarField.farFieldFromPowerPattern(PH(:),TH(:),P,freq,'linearY');
            else
                % Basic input error checking
                [Nang_x, Nf_x] = size(x);
                [Nang_y, Nf_y] = size(y);
                [Nang_E1, Nf_E1] = size(E1);
                [Nang_E2, Nf_E2] = size(E2);
                if nargin < 5 || isempty(E3)
                    E3 = zeros(size(E2));
                end
                [Nang_E3, Nf_E3] = size(E3);
                if ~isequal(Nang_x,Nang_y,Nang_E1,Nang_E2,Nang_E3)
                    error('All inputs must have the same number of rows');
                end
                if ~isequal(Nf_E1,Nf_E2,Nf_E3)
                    error('E1, E2, and E3 must have the same number of columns');
                end
                Nf = Nf_E1;
                if nargin < 6
                    freq = ones(1,Nf).*obj.freq;
                    %                 warning('freq not specified - using default for all columns');
                else
                    if isempty(freq)
                        freq = ones(1,Nf).*obj.freq;
                        %                     warning('freq not specified - using default for all columns');
                    elseif Nf ~= length(freq(1,:))
                        error('freq must have the same number of columns as E1 and E2');
                    end
                end
                if nargin < 7
                    Prad = ones(1,Nf).*obj.Prad;
                    %                 warning('Prad not specified - using default for all columns');
                else
                    if isempty(Prad)
                        Prad = ones(1,Nf).*obj.Prad;
                        %                     warning('Prad not specified - using default for all columns');
                    elseif Nf ~= length(Prad(1,:))
                        error('Prad must have the same number of columns as E1 and E2');
                    end
                end
                if nargin < 8
                    radEff = ones(1,Nf).*obj.radEff;
                    %                 warning('radEff not specified - using default for all columns');
                else
                    if isempty(radEff)
                        radEff = ones(1,Nf).*obj.radEff;
                        %                     warning('radEff not specified - using default for all columns');
                    elseif Nf ~= length(radEff(1,:))
                        error('radEff must have the same number of columns as E1 and E2');
                    end
                end
                if nargin < 9
                    coorType = obj.coorType;
                else
                    if isempty(coorType)
                        coorType = obj.coorType;
                    end
                end
                if nargin < 10
                    polType = obj.polType;
                else
                    if isempty(polType)
                        polType = obj.polType;
                    end
                end
                if nargin < 11
                    gridType = obj.gridType;
                else
                    if isempty(gridType)
                        gridType = obj.gridType;
                    end
                end
                if nargin < 12
                    freqUnit = obj.freqUnit;
                else
                    if isempty(freqUnit)
                        freqUnit = obj.freqUnit;
                    end
                end
                if nargin < 13
                    slant = obj.slant;
                else
                    if isempty(slant)
                        slant = obj.slant;
                    end
                end
                if Nf_x > 1 || Nf_y > 1
                    warning('Only using first column of th and ph since they must be equal for all frequencies anyway');
                    x = x(:,1);
                    y = y(:,1);
                end
                
                obj.x = x;
                obj.y = y;
                obj.E1 = E1;
                obj.E2 = E2;
                obj.E3 = E3;
                obj = setFreq(obj,freq,freqUnit);
                obj.Prad = Prad;
                obj.radEff = radEff;
                obj.coorType = coorType;
                obj.polType = polType;
                obj.gridType = gridType;
                obj.slant = slant;
                obj = setBase(obj);
            end
        end
        
        %% Setters
        function ph = get.ph(obj)
           [ph,~] = getPhThCurrent(obj);
        end
        
        function th = get.th(obj)
           [~,th] = getPhThCurrent(obj);
        end
        
        function Nf = get.Nf(obj)
            Nf = numel(obj.freq);
        end
        
        function Nx = get.Nx(obj)
            Nx = length(unique(obj.x));
        end
        
        function Ny = get.Ny(obj)
            Ny = length(unique(obj.y));
        end
        
        function Nang = get.Nang(obj)
            Nang = size(obj.x,1);
        end
        
        function freqHz = get.freqHz(obj)
            switch obj.freqUnit
                case 'Hz'
                    freqMult = 1;
                case 'kHz'
                    freqMult = 1e3;
                case 'MHz'
                    freqMult = 1e6;
                case 'GHz'
                    freqMult = 1e9;
                case 'THz'
                    freqMult = 1e12;
            end
            freqHz = obj.freq*freqMult;
        end
        
        function Directivity_dBi = get.Directivity_dBi(obj)
           Directivity_dBi = dB10(max(obj.getDirectivity())); 
        end
        
        function Gain_dB = get.Gain_dB(obj)
            Gain_dB = dB10(max(obj.getGain()));
        end
        
        function radEff_dB = get.radEff_dB(obj)
            radEff_dB = dB10(obj.radEff);
        end
        
        function xname = get.xname(obj)
           [xname,~] = setXYnames(obj);
        end
        
        function yname = get.yname(obj)
           [~,yname] = setXYnames(obj);
        end
        
        function E1name = get.E1name(obj)
           [E1name,~] = setEnames(obj);
        end
        
        function E2name = get.E2name(obj)
           [~,E2name] = setEnames(obj);
        end
        
        function xRangeType = get.xRangeType(obj)
            [xRangeType,~] = setRangeTypes(obj);
        end

        function yRangeType = get.yRangeType(obj)
            [~,yRangeType] = setRangeTypes(obj);
        end

        function obj = setFreq(obj,freq,freqUnit)
            if nargin > 1
                assert(numel(freq) == size(obj.E1,2),'Error, freq must be the same length as the number of columns in E1')
                obj.freq = freq;
            end
            if nargin > 2
                obj.freqUnit = freqUnit;
            end
        end
        
        %% Pattern getters
        function FFpattern = getFarFieldStruct(obj)
            % This returns the legacy structure format for testing with all
            % the tons of old code
            obj = obj.coor2spherical(true);
            FFpattern.th = repmat(obj.y,1,obj.Nf);
            FFpattern.ph = repmat(obj.x,1,obj.Nf);
            FFpattern.Eth = obj.E1;
            FFpattern.Eph = obj.E2;
            FFpattern.freq = obj.freqHz;
            FFpattern.Nth = obj.Ny;
            FFpattern.Nph = obj.Nx;
            FFpattern.Nf = obj.Nf;
            FFpattern.Prad = obj.Prad;
        end
        
        function [E1field, E2field, E3field] = getEfield(obj)
            % function [E1field, E2field, E3field] = getEfield(obj)
            % Returns the Efield matrices of size [Nang x Nf]
            % Efield = E*exp(-jkr)/r
            k = 2.*pi.*obj.freqHz./obj.c0;
            FFfact = exp(-1i.*k.*obj.r)./obj.r;
            E1field = bsxfun(@times,obj.E1,FFfact);
            E2field = bsxfun(@times,obj.E2,FFfact);
            E3field = bsxfun(@times,obj.E3,FFfact);
        end
        
        function [W] = getW(obj)
            % function [W] = getW(obj)
            % returns the radiation density in W/m2 [Nang x Nf]
            [E1f, E2f] = getEfield(obj);    % Can use any orthogonal pair
            W = 1./(2.*obj.eta0).*(abs(E1f).^2 + abs(E2f).^2);
        end
        
        function [U] = getU(obj)
            % function [U] = getU(obj)
            % returns the radiation intensity in W/unit solid angle [Nang x Nf]
            U = obj.r^2.*getW(obj);
        end
        
        function [D] = getDirectivity(obj)
            % function [D] = getDirectivity(obj)
            % returns the directivity (linear) in D [Nang x Nf]
            D = 4.*pi.*bsxfun(@times,getU(obj),1./obj.Prad);
        end
        
        function [G] = getGain(obj)
            % function [G] = getGain(obj)
            % returns the gain (linear) in G [Nang x Nf]
            G = bsxfun(@times,getDirectivity(obj),obj.radEff);
        end
        
        function [AR] = getAxialRatio(obj)
            % function [AR] = getAxialRatio(obj)
            % returns the Axial Ratio (linear) [Nang x Nf]
            AR = sqrt((abs(obj.E1).^2 + abs(obj.E2).^2 + abs(obj.E1.^2 + obj.E2.^2))./(abs(obj.E1).^2 + abs(obj.E2).^2 - abs(obj.E1.^2 + obj.E2.^2)));
        end
        
        function [ARinv] = getAxialRatioInv(obj)
            % function [ARinv] = getAxialRatioInv(obj)
            % returns the inverted Axial Ratio in ARinv [Nang x Nf]
            ARinv = sqrt((abs(obj.E1).^2 + abs(obj.E2).^2 - abs(obj.E1.^2 + obj.E2.^2))./(abs(obj.E1).^2 + abs(obj.E2).^2 + abs(obj.E1.^2 + obj.E2.^2)));
        end
        
        function [Xpol] = getCO_XP(obj)
            % function [Xpol] = getCO_XP(obj)
            % returns the CO/XP ratio (linear) [Nang x Nf]
            Xpol = (abs(obj.E2)./abs(obj.E1)).^2;
        end
        
        function [Xpol] = getXP_CO(obj)
            % function [Xpol] = getXP_CO(obj)
            % returns the XP/CO ratio (linear) [Nang x Nf]
            Xpol = (abs(obj.E1)./abs(obj.E2)).^2;
        end
        
        %% Grid getters
        function [u, v, w] = getDirCos(obj)
            handle2DirCos = str2func([obj.gridTypeBase,'2DirCos']);
            [u,v,w] = handle2DirCos(obj.xBase,obj.yBase);
        end
        
        function [ph, th] = getPhTh(obj)
            switch obj.gridTypeBase
                case 'PhTh'
                    ph = obj.xBase;
                    th = obj.yBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [ph,th] = DirCos2PhTh(u,v,w);
            end
        end
        
        function [az, el] = getAzEl(obj)
            switch obj.gridTypeBase
                case 'AzEl'
                    el = obj.yBase;
                    az = obj.xBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [az,el] = DirCos2AzEl(u,v,w);
            end
        end
        
        function [ep, al] = getElAz(obj)
            switch obj.gridTypeBase
                case 'ElAz'
                    ep = obj.xBase;
                    al = obj.yBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [ep,al] = DirCos2ElAz(u,v,w);
            end
        end
        
        function [Xg, Yg] = getTrueView(obj)
            switch obj.gridTypeBase
                case 'TrueView'
                    Xg = obj.xBase;
                    Yg = obj.xBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [Xg,Yg] = DirCos2TrueView(u,v,w);
            end
        end
        
        function [asinu, asinv] = getArcSin(obj)
            switch obj.gridTypeBase
                case 'ArcSin'
                    asinu = obj.xBase;
                    asinv = obj.yBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [asinu,asinv] = DirCos2ArcSin(u,v,w);
            end
        end
        
        %% Grid transformation setters
        function obj = grid2PhTh(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'PhTh')
                [obj.x,obj.y] = getPhTh(obj);
                obj.gridType = 'PhTh';
            end
        end
        
        function obj = grid2DirCos(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'DirCos')
                [obj.x,obj.y] = getDirCos(obj);
                obj.gridType = 'DirCos';
            end
        end
        
        function obj = grid2AzEl(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'AzEl')
                [obj.x,obj.y] = getAzEl(obj);
                obj.gridType = 'AzEl';
            end
        end
        
        function obj = grid2ElAz(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'ElAz')
                [obj.x,obj.y] = getElAz(obj);
                obj.gridType = 'ElAz';
            end
        end
        
        function obj = grid2TrueView(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'TrueView')
                [obj.x,obj.y] = getTrueView(obj);
                obj.gridType = 'TrueView';
            end
        end
        
        function obj = grid2ArcSin(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'ArcSin')
                [obj.x,obj.y] = getArcSin(obj);
                obj.gridType = 'ArcSin';
            end
        end
        
        %% Grid range shifters
        function obj = sortGrid(obj)
            obj = roundGrid(obj);
            [~,iSort] = sortrows([obj.x,obj.y],[1 2]);
            obj.x = obj.x(iSort);
            obj.y = obj.y(iSort);
            obj.E1 = obj.E1(iSort,:);
            obj.E2 = obj.E2(iSort,:);
            obj.E3 = obj.E3(iSort,:);
        end
        
        function obj = roundGrid(obj,nSigDig)
            % Round to some significant digits for sorting (some issues can
            % arise in deg2rad and rad2deg
            if nargin < 2
                nSigDig = obj.nSigDig;
            end
            xRound = round(obj.x*10^nSigDig)/10^nSigDig;
            yRound = round(obj.y*10^nSigDig)/10^nSigDig;
            obj.x = xRound;
            obj.y = yRound;
        end
        
        function obj = copyAndInsertXcut(obj1,xvalCopy,xvalPaste,tol)
            % Use this to copy an X cut into another position.  Typically
            % handy when some transformation does not include the closing
            % cut - that is the 0 and 360 or -180 and 180 cuts.  Can in
            % principle be used to do random stuff - so careful.
            
            if nargin < 4
                tol = mean(diff(unique(obj1.x)));
            end
            % Make a whole new object to initialise the base
            % correctly
            inInd = find(abs(obj1.x - xvalCopy) < tol);
            xNew = [obj1.x;xvalPaste.*ones(size(inInd))];
            yNew = [obj1.y;obj1.y(inInd)];
            E1New = [obj1.E1;obj1.E1(inInd,:)];
            E2New = [obj1.E2;obj1.E2(inInd,:)];
            E3New = [obj1.E3;obj1.E3(inInd,:)];
            obj = FarField(xNew,yNew,E1New,E2New,E3New,obj1.freq,...
                obj1.Prad.*2,obj1.radEff,obj1.coorType,obj1.polType,obj1.gridType,obj1.freqUnit,obj1.slant);
            obj = obj.sortGrid;
            obj = FarField(obj.x,obj.y,obj.E1,obj.E2,obj.E3,obj.freq,...
                obj.Prad.*2,obj.radEff,obj.coorType,obj.polType,obj.gridType,obj.freqUnit,obj.slant);
        end
        
        obj = setXrange(obj,type)
        obj = setYrange(obj,type)
        
        %% Coordinate system getters
        function [Eth, Eph, Er] = getEspherical(obj)
            [Ph,Th] = getPhTh(obj);
            TH = repmat(Th(:,1),1,obj.Nf);
            PH = repmat(Ph(:,1),1,obj.Nf);
            % Change to the Base values here....
            switch obj.coorTypeBase
                case 'spherical'
                    Eth = obj.E1Base;
                    Eph = obj.E2Base;
                case 'Ludwig1'
                    Eth = cos(TH).*cos(PH).*obj.E1Base + cos(TH).*sin(PH).*obj.E2Base - sin(TH).*obj.E3Base;
                    Eph = -sin(PH).*obj.E1Base + cos(PH).*obj.E2Base;
                case 'Ludwig2AE'
                    cosEl = sqrt(1 - sin(TH).^2.*sin(PH).^2);
                    Del = cos(PH).^2 + cos(TH).^2.*sin(PH).^2;
                    Eth = (cosEl./Del).*(cos(PH).*obj.E1Base + cos(TH).*sin(PH).*obj.E2Base);
                    Eph = (cosEl./Del).*(-cos(TH).*sin(PH).*obj.E1Base + cos(PH).*obj.E2Base);
                case 'Ludwig2EA'
                    cosAl = sqrt(1 - sin(TH).^2.*cos(PH).^2);
                    Del = cos(TH).^2.*cos(PH).^2 + sin(PH).^2;
                    Eth = (cosAl./Del).*(cos(TH).*cos(PH).*obj.E1Base + sin(PH).*obj.E2Base);
                    Eph = (cosAl./Del).*(-sin(PH).*obj.E1Base + cos(TH).*cos(PH).*obj.E2Base);
                case 'Ludwig3'
                    Del = 1;
                    Eth = (1./Del).*(cos(PH).*obj.E1Base + sin(PH).*obj.E2Base);
                    Eph = (1./Del).*(-sin(PH).*obj.E1Base + cos(PH).*obj.E2Base);
            end
            Er = zeros(size(Eth));
        end
        
        function [Ex, Ey, Ez] = getELudwig1(obj)
            switch obj.coorTypeBase
                case 'Ludwig1'
                    Ex = obj.E1Base;
                    Ey = obj.E2Base;
                    Ez = obj.E3Base;
                otherwise
                    [Eth, Eph, ~] = getEspherical(obj);
                    TH = repmat(obj.th(:,1),1,obj.Nf);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    % Assume farfield so no radial E-field
                    Ex = cos(TH).*cos(PH).*Eth - sin(PH).*Eph;
                    Ey = cos(TH).*sin(PH).*Eth + cos(PH).*Eph;
                    Ez = zeros(size(Ex));   % Strict definition in the Ludwig paper
            end
        end
        
        function [Eaz, Eel, E3] = getELudwig2AE(obj)
            switch obj.coorTypeBase
                case 'Ludwig2AE'
                    Eaz = obj.E1Base;
                    Eel = obj.E2Base;
                    E3 = obj.E3Base;
                otherwise
                    [Eth, Eph] = getEspherical(obj);
                    TH = repmat(obj.th(:,1),1,obj.Nf);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    cosEl = sqrt(1 - sin(TH).^2.*sin(PH).^2);
                    Eaz = (1./cosEl).*(cos(PH).*Eth - cos(TH).*sin(PH).*Eph);
                    Eel = (1./cosEl).*(cos(TH).*sin(PH).*Eth + cos(PH).*Eph);
                    E3 = zeros(size(Eaz));
                    % Sort out singularities poles - ToDo
                    phPoles = deg2rad([-270,-90,90,270].');
                    poleMat = [ones(4,1).*deg2rad(90),phPoles]; % [th=90,ph]
                    [~,iPole] = ismember(poleMat,[obj.th,obj.ph],'rows');
                    iPole = iPole(iPole>0);
                    [Eaz(iPole,:),Eel(iPole,:),E3(iPole,:)] = deal(0);
                    %                     Eaz(iPole,:) = Eth(iPole,:) + Eph(iPole,:);
            end
        end
        
        function [Eal, Eep, E3] = getELudwig2EA(obj)
            switch obj.coorTypeBase
                case 'Ludwig2EA'
                    Eal = obj.E1Base;
                    Eep = obj.E2Base;
                    E3 = obj.E3Base;
                otherwise
                    [Eth, Eph] = getEspherical(obj);
                    TH = repmat(obj.th(:,1),1,obj.Nf);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    cosAl = sqrt(1 - sin(TH).^2.*cos(PH).^2);
                    Eal = (1./cosAl).*(cos(TH).*cos(PH).*Eth - sin(PH).*Eph);
                    Eep = (1./cosAl).*(sin(PH).*Eth + cos(TH).*cos(PH).*Eph);
                    E3 = zeros(size(Eal));
                    % Sort out singularities poles - ToDo
                    phPoles = deg2rad([-360,-180,0,180,360].');
                    poleMat = [ones(5,1).*deg2rad(90),phPoles]; % [th=90,ph]
                    [~,iPole] = ismember(poleMat,[obj.th,obj.ph],'rows');
                    iPole = iPole(iPole>0);
                    [Eal(iPole,:),Eep(iPole,:),E3(iPole,:)] = deal(0);
            end
        end
        
        function [Eh, Ev, E3] = getELudwig3(obj)
            switch obj.coorTypeBase
                case 'Ludwig3'
                    Eh = obj.E1Base;
                    Ev = obj.E2Base;
                    E3 = obj.E3Base;
                otherwise
                    [Eth, Eph] = getEspherical(obj);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    Eh = cos(PH).*Eth - sin(PH).*Eph;
                    Ev = sin(PH).*Eth + cos(PH).*Eph;
                    E3 = zeros(size(Eh));
            end
        end
        
        
        %% Coordinate system transformation methods
        function obj = coor2spherical(obj,setStdGrid)
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'spherical')
                [obj.E1,obj.E2,obj.E3] = getEspherical(obj);
                obj.E3 = zeros(size(obj.E1));
                obj.coorType = 'spherical';
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        function obj = coor2Ludwig1(obj,setStdGrid)
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig1')
                [obj.E1,obj.E2,obj.E3] = getELudwig1(obj);
                obj.coorType = 'Ludwig1';
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        function obj = coor2Ludwig2AE(obj,setStdGrid)
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig2AE')
                [obj.E1,obj.E2,obj.E3] = getELudwig2AE(obj);
                obj.coorType = 'Ludwig2AE';
            end
            if setStdGrid
                obj = obj.grid2AzEl;
            end
        end
        
        function obj = coor2Ludwig2EA(obj,setStdGrid)
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig2EA')
                [obj.E1,obj.E2,obj.E3] = getELudwig2EA(obj);
                obj.coorType = 'Ludwig2EA';
            end
            if setStdGrid
                obj = obj.grid2ElAz;
            end
        end
        
        function obj = coor2Ludwig3(obj,setStdGrid)
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig3')
                [obj.E1,obj.E2,obj.E3] = getELudwig3(obj);
                obj.coorType = 'Ludwig3';
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        %% Polarization type getters
        function [E1lin, E2lin, E3lin] = getElin(obj)
            
            % Start at the base, transform to correct coordinate system
            coorTypeIn = obj.coorType;
            coorTypeH = str2func(['coor2',coorTypeIn]);
            obj1 = obj.reset2Base;
            obj1 = coorTypeH(obj1,false);
            switch obj.polTypeBase % Should be the same as the transformed object - can use obj or obj1
                case 'linear'
                    E1lin = obj1.E1;
                    E2lin = obj1.E2;
                case 'circular'
                    Del = 2*1i;
                    E1lin = sqrt(2)./Del.*(1i.*obj1.E1 + 1i.*obj1.E2);
                    E2lin = sqrt(2)./Del.*(-obj1.E1 + obj1.E2);
                case 'slant'
                    PSI = ones(size(obj1.E1)).*obj.slant; % Slant of input object
                    Del = 1;
                    E1lin = 1./Del.*(cos(PSI).*obj1.E1 + sin(PSI).*obj1.E2);
                    E2lin = 1./Del.*(-sin(PSI).*obj1.E1 + cos(PSI).*obj1.E2);
            end
            E3lin = zeros(size(E1lin));
        end
        
        function [Elh,Erh,E3circ] = getEcircular(obj)
            switch obj.polType
                case 'circular'
                    Elh = obj.E1;
                    Erh = obj.E2;
                    E3circ = obj.E3;
                otherwise
                    [E1lin, E2lin] = getElin(obj);
                    Elh = 1/sqrt(2).*(E1lin - 1i.*E2lin);
                    Erh = 1/sqrt(2).*(E1lin + 1i.*E2lin);
                    E3circ = zeros(size(Elh));
            end
        end
        
        function [Exp,Eco,E3slant] = getEslant(obj)
            switch obj.polType
                case 'slant'
                    Exp = obj.E1;
                    Eco = obj.E2;
                    E3slant = obj.E3;
                otherwise
                    [E1lin, E2lin] = getElin(obj);
                    PSI = ones(size(obj.E1)).*obj.slant;
                    Exp = cos(PSI).*E1lin - sin(PSI).*E2lin;
                    Eco = sin(PSI).*E1lin + cos(PSI).*E2lin;
                    E3slant = zeros(size(Exp));
            end
        end
        
        
        %% Polarisation type transformation methods
        function obj = pol2linear(obj)
            if ~strcmp(obj.polType,'linear')
                [obj.E1, obj.E2, obj.E3] = getElin(obj);
                obj.polType = 'linear';
            end
        end
        
        function obj = pol2circular(obj)
            if ~strcmp(obj.polType,'circular')
                [obj.E1,obj.E2,obj.E3] = getEcircular(obj);
                obj.polType = 'circular';
            end
        end
        
        function obj = pol2slant(obj)
            if ~strcmp(obj.polType,'slant')
                [obj.E1,obj.E2,obj.E3] = getEslant(obj);
                obj.polType = 'slant';
            end
        end
        
        %% Format transformation
        function obj = transformTypes(obj, obj1)
            % Function to transform the format of obj to that of obj1 -
            % that is the grid,coor, and pol Types of obj goes to those of
            % obj1.
            objGridType = obj1.gridType;
            objCoorType = obj1.coorType;
            objPolType = obj1.polType;
            handleGridType = str2func(['grid2',objGridType]);
            handleCoorType = str2func(['coor2',objCoorType]);
            handlePolType = str2func(['pol2',objPolType]);
            obj = handleGridType(obj);
            obj = handleCoorType(obj,0);
            obj = handlePolType(obj);
        end

        %% Base grid functions
        function obj = reset2Base(obj)
            % Hard reset to the base format
            obj.x = obj.xBase;
            obj.y = obj.yBase;
            obj.gridType = obj.gridTypeBase;
            obj.E1 = obj.E1Base;
            obj.E2 = obj.E2Base;
            obj.E3 = obj.E3Base;
            obj.coorType = obj.coorTypeBase;
            obj.polType = obj.polTypeBase;
%             obj = setRangeTypes(obj);
        end
        
        function obj = grid2Base(obj)
            % Evaluate the current object (pol and coor) on the base grid
            coorTypeIn = obj.coorType;
            polTypeIn = obj.polType;
            coorTypeH = str2func(['coor2',coorTypeIn]);
            polTypeH = str2func(['pol2',polTypeIn]);
            obj = obj.reset2Base;
            % Keep the current coorType and polType
            obj = coorTypeH(obj,false);
            obj = polTypeH(obj);
%             obj = setRangeTypes(obj);
        end
        
        function obj = currentForm2Base(obj1,stepDeg,xylimsDeg,hemisphere)
            % Sets the base to the current format. Resamples the field on a
            % regular plaid grid, and makes this the new base grid. This is
            % typically then not where actual samples where, but instead
            % interpolated values. 
            % Format of xylimsDeg is [xmin xmax; ymin ymax]
            
            % Due to the interpolation only operating on the y in 180 range
            % from the Direction Cosine transforms...
            assert(~isequal(obj1.yRangeType,'360'),'yRangeType cannot be 360 for the base representation')
            
            % Set defaults
            if strcmp(obj1.gridTypeBase,'DirCos') || strcmp(obj1.gridType,'ArcSin')
                stepX = asin(min(abs(diff(unique(obj1.xBase)))));
                stepY = asin(min(abs(diff(unique(obj1.yBase)))));
                xmin = min(asin(obj1.x));
                xmax = max(asin(obj1.x));
                ymin = min(asin(obj1.y));
                ymax = max(asin(obj1.y));
            else
                % Sort out rounding errors for degrees
                stepX = deg2rad(round(rad2deg(min(abs(diff(unique(obj1.xBase)))))*10^obj1.nSigDig)/10^obj1.nSigDig);
                stepY = deg2rad(round(rad2deg(min(abs(diff(unique(obj1.yBase)))))*10^obj1.nSigDig)/10^obj1.nSigDig);
                xmin = deg2rad(round(rad2deg(min(obj1.x))*10^obj1.nSigDig)./10^obj1.nSigDig);
                xmax = deg2rad(round(rad2deg(max(obj1.x))*10^obj1.nSigDig)./10^obj1.nSigDig);
                ymin = deg2rad(round(rad2deg(min(obj1.y))*10^obj1.nSigDig)./10^obj1.nSigDig);
                ymax = deg2rad(round(rad2deg(max(obj1.y))*10^obj1.nSigDig)./10^obj1.nSigDig);
            end
            
            hem = 'top';
            
            % Overwrite defaults
            if nargin >= 2
                if numel(stepDeg) == 1
                    [stepX,stepY] = deal(deg2rad(stepDeg));
                elseif numel(stepDeg) == 2
                    stepX = deg2rad(stepDeg(1));
                    stepY = deg2rad(stepDeg(2));
                end
                if nargin >= 3
                    xylim = deg2rad(xylimsDeg);
                    xmin = xylim(1,1);
                    xmax = xylim(1,2);
                    ymin = xylim(2,1);
                    ymax = xylim(2,2);
                    if nargin == 4
                        hem = hemisphere;
                    end
                end
            end
            
            if strcmp(obj1.gridType,'DirCos') || strcmp(obj1.gridType,'ArcSin')
                stepX = sin(stepX);
                stepY = sin(stepY);
                xmin = sin(xmin);
                xmax = sin(xmax);
                ymin = sin(ymin);
                ymax = sin(ymax);
            end
            
            % Build the new grid
            Nxi = round((xmax - xmin)/stepX) + 1;
            Nyi = round((ymax - ymin)/stepY) + 1;
            xivect = linspace(xmin,xmax,Nxi);
            yivect = linspace(ymin,ymax,Nyi);
            [Xi,Yi] = meshgrid(xivect,yivect);
            xi = Xi(:);
            yi = Yi(:);

            % Interpolate the fields
            [E1grid,E2grid,E3grid] = deal(zeros(Nxi*Nyi,obj1.Nf));
            for ff = 1:obj1.Nf
                E1grid(:,ff) = interpolateGrid(obj1,'E1',xi,yi,ff,hem);
                E2grid(:,ff) = interpolateGrid(obj1,'E2',xi,yi,ff,hem);
                E3grid(:,ff) = interpolateGrid(obj1,'E3',xi,yi,ff,hem);
            end
            % Remove the extra phase introduced by the interpolateGrid
            % function - this just keeps the real/imag and phase field
            % consistant with the plotting
            k = 2.*pi.*obj1.freqHz./obj1.c0;
            FFfact = exp(1i.*k.*obj1.r)./obj1.r;
            E1grid = bsxfun(@times,E1grid,FFfact);
            E2grid = bsxfun(@times,E2grid,FFfact);
            E3grid = bsxfun(@times,E3grid,FFfact);
            
            % Populate the new farField object
            obj = obj1;
            obj.x = xi(:);
            obj.y = yi(:);
            obj.E1 = E1grid;
            obj.E2 = E2grid;
            obj.E3 = E3grid;
            obj = setBase(obj);
            % Update the current form to the base form
            obj = reset2Base(obj);
        end
        
        
        %% Plotting methods
        plot(obj,varargin)
        plotJones(obj1,obj2,varargin)  % Much to do here still...
        plotPrincipleCuts(obj,varargin)
        
        function plotGrid(obj,markerStyle)
            if nargin < 2
                markerStyle = 'k.';
            end
            switch obj.gridType
                case {'DirCos'}
                    xplot = obj.x;
                    yplot = obj.y;
                    xtext = [obj.xname, ' = sin(\theta)cos(\phi)'];
                    ytext = [obj.yname, ' = sin(\theta)sin(\phi)'];
                    
                otherwise
                    xplot = rad2deg(obj.x);
                    yplot = rad2deg(obj.y);
                    xtext = [obj.xname,' (deg)'];
                    ytext = [obj.yname,' (deg)'];
            end
            plot(xplot,yplot,markerStyle)
            xlabel(xtext)
            ylabel(ytext)
            axis equal
            grid on
            xlim([min(xplot),max(xplot)])
            ylim([min(yplot),max(yplot)])
        end
        
        function plotGridBase(obj)
            obj = obj.reset2Base;
            obj.plotGrid;
        end
        
        %% Interpolation methods
        [Z] = interpolateGrid(obj,output,xi,yi,varargin)
        
        %% Phase centre/shifts/rotations of the field
        [Z, Delta, delta0, eta_pd] = phaseCentreKildal(FF,pol,th_M)

        function obj = rotate(obj1,rotHandle,rotAng)
            % General rotation function for FarField objects
            % rotHandle is the function handle for the type of rotation:
            %   rotx3Dsph, roty3Dsph, rotz3Dsph, rotGRASPsph, rotEulersph
            % rotAng is the associated angle in rad. Scalar for rotations
            % around an axis, and [3x1] for GRASP or Euler rotations
            
            % THis will probably depent on the pattern type which one is
            % best.  Only have spherical and Ludwig 3 implemented for now,
            % so hard-coded.
            baseCoorType = 'Ludwig3';
            coorHandle = str2func(['coor2',baseCoorType]);

            % Test if the rotation function handle has the trailing 'sph'
            handleStr = func2str(rotHandle);
            if ~strcmp(handleStr(end-2:end),'sph')
                handleStr = [handleStr,'sph'];  % Add it if not - some user errors fixed at least!
                rotHandle = str2func(handleStr);
            end
            % Transform to sensible grid and coordinate system for rotation
            FFsph = obj1.grid2PhTh;  % Always work in the PhTh coordinate system
            FFsph = FFsph.setXrange('sym'); % Always work in symmetrical xRange
            FFsph = coorHandle(FFsph,false);
            
%             % Force all the th = 0|180 fields to be identical - fixes pole
%             % interpolation problems
%             tol = 10^(-obj1.nSigDig);
%             i_0_0 = find(abs(FFsph.th) < tol & abs(FFsph.ph) < tol);
%             i_180_0 = find(abs(FFsph.th - pi) < tol & abs(FFsph.ph) < tol);
%             if ~isempty(i_0_0)
%                 E1_0_0 = obj1.E1(i_0_0(1),:);
%                 E2_0_0 = obj1.E2(i_0_0(1),:);
%                 i_0 = find(abs(FFsph.th - 0)<tol);
%                 FFsph.E1(i_0,:) = repmat(E1_0_0,length(i_0),1);
%                 FFsph.E2(i_0,:) = repmat(E2_0_0,length(i_0),1);
%             end
%             if ~isempty(i_180_0)
%                 E1_180_0 = obj1.E1(i_180_0(1),:);
%                 E2_180_0 = obj1.E2(i_180_0(1),:);
%                 i_180 = find(abs(FFsph.th - pi)<tol);
%                 FFsph.E1(i_180,:) = repmat(E1_180_0,length(i_180),1);
%                 FFsph.E2(i_180,:) = repmat(E2_180_0,length(i_180),1);
%             end

            % Get the grid step sizes from the original
            stepx = (max(FFsph.x) - min(FFsph.x))./(FFsph.Nx-1);
            stepy = (max(FFsph.y) - min(FFsph.y))./(FFsph.Ny-1);
            stepDeg = rad2deg([stepx,stepy]);
            xmin = min(FFsph.x);
            xmax = max(FFsph.x);
            ymin = min(FFsph.y);
            ymax = max(FFsph.y);
            % Perform the rotation of the grid
            phIn = FFsph.x.';
            thIn = FFsph.y.';
            sphAngIn = [phIn;thIn];
            sphAngRot = rotHandle(sphAngIn,rotAng);
            phOut = sphAngRot(1,:).';
            thOut = sphAngRot(2,:).';
            FFsph.x = phOut;
            FFsph.y = thOut;

            % Perform the rotation of the field vectors
            % Coordinate systems
            C0 = coordinateSystem;
            rotHandleStr = func2str(rotHandle);
            rotHandleCoor = str2func(rotHandleStr(1:end-3));
            Crot = rotHandleCoor(C0,rotAng);
            
            % Vector origin points before rotation
            [OIx,OIy,OIz] = PhTh2DirCos(phIn,thIn);
            origin_In = pnt3D(OIx(:).',OIy(:).',OIz(:).');
            origin_out = origin_In.changeBase(C0,Crot);

            % local unit vector directions before and after rotation
            switch baseCoorType
                case 'spherical'
                    [xHatIn,yHatIn,~] = unitVectorsSpherical(thIn,phIn);
                    [xHatOut,yHatOut,~] = unitVectorsSpherical(thOut,phOut);
                    [E1sign,E2sign] = deal(1,1);
                case 'Ludwig3'
                    [xHatIn,yHatIn] = unitVectorsDirCos(thIn,phIn);
                    [xHatOut,yHatOut] = unitVectorsDirCos(thOut,phOut);
                    [E1sign,E2sign] = deal(-1,-1);
            end
            % Vector tip points before rotation
            xTipIn = origin_In + xHatIn;
            yTipIn = origin_In + yHatIn;
            % Rotate all the points
            xTipOut = xTipIn.changeBase(C0,Crot);
            yTipOut = yTipIn.changeBase(C0,Crot);
            xOut = xTipOut - origin_out;
            yOut = yTipOut - origin_out;
            % Project onto the local unit vectors
            xx = dot(xOut.pointMatrix,xHatOut.pointMatrix);
            xy = dot(xOut.pointMatrix,yHatOut.pointMatrix);
            yx = dot(yOut.pointMatrix,xHatOut.pointMatrix);
            yy = dot(yOut.pointMatrix,yHatOut.pointMatrix);
            E1rot = FFsph.E1.*repmat(xx(:),1,FFsph.Nf) + FFsph.E2.*repmat(yx(:),1,FFsph.Nf);
            E2rot = FFsph.E1.*repmat(xy(:),1,FFsph.Nf) + FFsph.E2.*repmat(yy(:),1,FFsph.Nf);
            FFsph.E1 = E1sign.*E1rot;
            FFsph.E2 = E2sign.*E2rot;

            % Set the baseGrid of the rotated object.  This is required
            % since all transformations operate from the base grid
            FFsph = FFsph.sortGrid;
            FFsph = FFsph.setBase;
            FFsph = FFsph.currentForm2Base(stepDeg,rad2deg([xmin,xmax;ymin,ymax]));
            % Reset the grid and coordinate system, and reset the base back
            % in the original format
            obj = transformTypes(FFsph, obj1);
            obj = obj.setXrange(obj1.xRangeType);
            obj = obj.currentForm2Base();
        end
        
        function obj = shift(obj,shiftVect)
           % Shifts the FarField by a distance specified in the pnt3D input
           % shiftVect (only uses the first entry)
           % shiftVect can also be a vector of length 3 with elements
           % [delX,delY,delZ] in m
           if nargin < 2, shiftVect = pnt3D; end
           if ~isa(shiftVect,'pnt3D')
               assert(numel(shiftVect)==3,'Expect a vector of length 3 for the shifVect, if not a pnt3D');
               shiftVect = pnt3D(shiftVect(1),shiftVect(2),shiftVect(3));
           end
           assert(numel(shiftVect.x)==1,'Only one point allowed when shifting a FarField');
           
           lam = obj.c0./obj.freqHz;
           k = 2.*pi./lam;
           
           r_hat = [sin(obj.th).*cos(obj.ph), sin(obj.th).*sin(obj.ph), cos(obj.th)];
           rmat = repmat(shiftVect.pointMatrix.',obj.Nang,1);
           rdotr = dot(rmat,r_hat,2);
           phase = exp(1i.*bsxfun(@times,k,rdotr));
           obj.E1 = obj.E1.*phase;
           obj.E2 = obj.E2.*phase;
           obj.E3 = obj.E3.*phase;
           obj = obj.setBase;
        end
        
        
        %% Maths
        function obj = plus(obj1,obj2)
            obj1base = reset2Base(obj1);
            obj2base = reset2Base(obj2);
            
            if isGridEqual(obj1base,obj2base) && typesAreEqual(obj1base,obj2base)
                obj = obj1base;
                obj.E1 = obj1base.E1 + obj2base.E1;
                obj.E2 = obj1base.E2 + obj2base.E2;
                obj.E3 = obj1base.E3 + obj2base.E3;
                obj.Prad = obj1base.Prad + obj2base.Prad;
                Pt = obj1base.Prad./obj1base.radEff + obj2base.Prad./obj2base.radEff;
                obj.radEff = obj.Prad./Pt;
                obj = setBase(obj);
            else
                error('Can only add FarFields with equal base grids')
            end
            
            if typesAreEqual(obj1,obj2)
                obj = transformTypes(obj, obj1); 
            end
        end
        
        function obj = minus(obj1,obj2)
            obj1base = reset2Base(obj1);
            obj2base = reset2Base(obj2);
            
            if isGridEqual(obj1base,obj2base) && typesAreEqual(obj1base,obj2base)
                obj = obj1base;
                obj.E1 = obj1base.E1 - obj2base.E1;
                obj.E2 = obj1base.E2 - obj2base.E2;
                obj.E3 = obj1base.E3 - obj2base.E3;
                obj.Prad = obj1base.Prad + obj2base.Prad;
                Pt = obj1base.Prad./obj1base.radEff + obj2base.Prad./obj2base.radEff;
                obj.radEff = obj.Prad./Pt;
                obj = setBase(obj);
            else
                error('Can only subtract FarFields with equal base grids')
            end
            
             if typesAreEqual(obj1,obj2)
                obj = transformTypes(obj, obj1);   
            end
        end
        
        function obj = times(obj1,obj2)
            % Don't operate in BaseGrid - straight on the actual values.
            % This is often used with conj, which operates in the current
            % grid.
            if isGridEqual(obj1,obj2) && typesAreEqual(obj1,obj2)
                obj = obj1;
                obj.E1 = obj1.E1.*obj2.E1;
                obj.E2 = obj1.E2.*obj2.E2;
                obj.E3 = obj1.E3.*obj2.E3;
                obj.Prad = obj.pradInt;
                obj.radEff = ones(size(obj.Prad));
                obj = setBase(obj);
            else
                error('Can only multiply FarFields with equal grids')
            end
        end
        
        function obj = conj(obj1)
            obj = obj1;
            obj.E1 = conj(obj1.E1);
            obj.E2 = conj(obj1.E2);
            obj.E3 = conj(obj1.E3);  
        end
        
        function obj = abs(obj1)
            obj = obj1;
            obj.E1 = abs(obj1.E1);
            obj.E2 = abs(obj1.E2);
            obj.E3 = abs(obj1.E3);  
        end
        
        function obj = scale(obj1,scaleFactor)
            % Scale the FarField object E-fields by the scaleFactor
            obj = obj1;
            obj.E1 = obj1.E1.*scaleFactor;
            obj.E2 = obj1.E2.*scaleFactor;
            obj.E3 = obj1.E3.*scaleFactor;
            obj.Prad = obj1.Prad.*(abs(scaleFactor).^2);
        end
        
        function [normE] = norm(obj,Ntype)
           % Calculate the vector norm of the three E-field components - overloads the MATLAB norm function
           % This is used for error checking mostly - when comparing
           % different fields for instance...
           if nargin == 1
               Ntype = 2;
           end
           nE1 = norm(obj.E1,Ntype);
           nE2 = norm(obj.E2,Ntype);
           nE3 = norm(obj.E3,Ntype);
           normE = [nE1,nE2,nE3];
        end
        
        function [rmsE1,rmsE2,rmsE3] = rms(obj,DIM)
           % Calculate the rms of the three E-field components - overloads the MATLAB rms function
           % This is used for error checking mostly - when comparing
           % different fields for instance...
           if nargin < 2
               DIM = 1;
           end
           rmsE1 = rms(obj.E1,DIM);
           rmsE2 = rms(obj.E2,DIM);
           rmsE3 = rms(obj.E3,DIM);
        end
        
        function T = convPower(obj1,obj2)
            % Convolve the power patterns, over the full sphere, of two
            % FarField objects. Typically used for antenna temperature
            % calculations
            obj1 = reset2Base(obj1);
            obj2 = reset2Base(obj2);
            
            if isGridEqual(obj1,obj2)
                P = obj1.getU.*obj2.getU;
                FF_T = FarField.farFieldFromPowerPattern(obj1.phBase,obj1.thBase,P,obj1.freq);
                T = FF_T.pradInt;
            else
                error('Can only convolve FarFields with equal base grids')
            end
        end
        
        %% Field normalization
        function P = pradInt(obj)
            % Returns the total power in the field integrated over the
            % full available grid
            obj = reset2Base(obj);
            symFact = 2^(sum(abs([obj.symXY,obj.symXZ,obj.symYZ])));
            assert(obj.isGridUniform,'Must have a plaid, monotonic, uniform grid for power calculation through integration');
            switch obj.gridType
                case 'PhTh'
                    PH = reshape(obj.x,obj.Ny,obj.Nx);
                    TH = reshape(obj.y,obj.Ny,obj.Nx);
                    U = obj.getU;
                    P = zeros(1,obj.Nf);
                    for ff = 1:obj.Nf
                        if ~obj.symmetryBOR1
                            integrand = reshape(U(:,ff),obj.Ny,obj.Nx).*sin(TH);
                            P(ff) = integral2D(PH,TH,integrand);
                        else
                            Nth = obj.Ny;
                            th_vect = obj.y(1:Nth);
                            integrand = (U(1:Nth,ff) + U(Nth+1:end,ff)).*sin(th_vect);
                            P(ff) = pi*integral1D(th_vect,integrand);
                            symFact = 1;    % Just to be sure...
                        end
                    end
                otherwise
                    error(['pradInt not implemented for gridType = ',obj.gridType,', only for PhTh grids'])
            end
            P = P.*symFact;
        end
        
        function obj = setPower(obj1,powerWatt)
            % Normalizes the FarField object obj1 to have the a total
            % radiated power specified in powerWatt (in Watts, of course)
            % powerWatt can be a vector of length obj.Nf or scalar.
            % The field need not be specified over the full sphere - the
            % total intercepted power in the specified sector will be set
            % to powerWatt.  Default value of 4*pi W will be used for one
            % argument.
            % The grid should be the standard plaid, montonic, uniform grid.
            obj1 = reset2Base(obj1);
            if nargin == 1
                powerWatt = 4*pi;
            end
            if length(powerWatt) == 1
                powerWatt = repmat(powerWatt,1,obj1.Nf);
            elseif length(powerWatt) ~= obj1.Nf
                error('powerWatt should be scalar or of length obj1.Nf');
            end
            P = obj1.Prad;
            Cn = powerWatt./(P);
            obj = obj1;
            obj = scale(obj,sqrt(Cn));
        end
        
        %% Frequency modifications
        function obj = getFi(obj1,freqIndex)
            % Returns an object only containing the results in freqIndex
            obj = obj1;
            obj.E1 = obj1.E1(:,freqIndex);
            obj.E2 = obj1.E2(:,freqIndex);
            obj.E3 = obj1.E3(:,freqIndex);
            obj.freq = obj1.freq(freqIndex);
            obj.Prad = obj1.Prad(freqIndex);
            obj.radEff = obj1.radEff(freqIndex);
            obj = setBase(obj);
        end
        
        %% Symmetry handlers
        function obj = setSymmetryXZ(obj,symmetryType)
            % Test if the input range is valid
            tol = 10^(-obj.nSigDig);
            % Easy to check in TrueView
            obj1 = obj.grid2TrueView;
            assert(all(sign(obj1.y+tol) > 0) || all(sign(obj1.y-tol) < 0),'Invalid range for XZ symmetry')
            obj.symmetryXZ = symmetryType;
            switch symmetryType
                case 'none'
                    obj.symXZ = 0;
                case 'electric'
                    obj.symXZ = -1;
                case 'magnetic'
                    obj.symXZ = 1;
                otherwise
                    error(['Unknown symmetry setting: ',symmetryType])
            end
        end
        
        function obj = setSymmetryYZ(obj,symmetryType)
            % Test if the input range is valid
            tol = 10^(-obj.nSigDig);
            % Easy to check in TrueView
            obj1 = obj.grid2TrueView;
            assert(all(sign(obj1.x+tol) > 0) || all(sign(obj1.x-tol) < 0),'Invalid range for YZ symmetry')
            obj.symmetryYZ = symmetryType;
            switch symmetryType
                case 'none'
                    obj.symYZ = 0;
                case 'electric'
                    obj.symYZ = -1;
                case 'magnetic'
                    obj.symYZ = 1;
                otherwise
                    error(['Unknown symmetry setting: ',symmetryType])
            end
        end
        
        function obj = setSymmetryXY(obj,symmetryType)
           warning('function: setSymmetryXY not implemented yet - unchanged object returned'); 
        end
        
        function obj = mirrorSymmetricPattern(obj1)
            % Returns the full pattern mirrored according to the symmetry
            % definitions
            
            if ~obj1.symXZ && ~obj1.symYZ && ~obj1.symXZ
                obj = obj1;
            else
                gridTypeIn = obj1.gridType;
                coorTypeIn = obj1.coorType;
                if strcmp(obj1.gridTypeBase,'DirCos') || strcmp(obj1.gridType,'ArcSin')
                    stepX = asin(min(abs(diff(unique(obj1.xBase)))));
                    stepY = asin(min(abs(diff(unique(obj1.yBase)))));
                else
                    % Sort out rounding errors for degrees
                    stepX = deg2rad(round(rad2deg(min(abs(diff(unique(obj1.xBase)))))*10^obj1.nSigDig)/10^obj1.nSigDig);
                    stepY = deg2rad(round(rad2deg(min(abs(diff(unique(obj1.yBase)))))*10^obj1.nSigDig)/10^obj1.nSigDig);
                end
                
                gridHandle = str2func(['grid2',gridTypeIn]);
                coorHandle = str2func(['coor2',coorTypeIn]);
                obj1 = obj1.grid2TrueView;
                obj1 = obj1.coor2Ludwig3(false);   % Always work in H/V for symmetry calculations...
                
                % Initialise
                XIn = [obj1.x];
                YIn = [obj1.y];
                E1In = [obj1.E1];
                E2In = [obj1.E2];
                E3In = [obj1.E3];
                if obj1.symXZ
                    XIn = [XIn;XIn];
                    YIn = [YIn;-YIn];
                    E1In = [E1In;obj1.symXZ.*E1In]; % Mirror according to symmetry
                    E2In = [E2In;-obj1.symXZ.*E2In];  % Mirror according to symmetry
                    E3In = [E3In;E3In];  % Do nothing for FarFields...
                end
                if obj1.symYZ
                    XIn = [XIn;-XIn];
                    YIn = [YIn;YIn];
                    E1In = [E1In;-obj1.symYZ.*E1In];  % Mirror according to symmetry
                    E2In = [E2In;obj1.symYZ.*E2In]; % Mirror according to symmetry
                    E3In = [E3In;E3In];  % Do nothing for FarFields...
                end
                % Object for grid/coor transformation
                objD = FarField(XIn,YIn,E1In,E2In,E3In,...
                    obj1.freq,obj1.Prad.*2,obj1.radEff,obj1.coorType,obj1.polType,obj1.gridType,obj1.freqUnit,obj1.slant);
                obj = coorHandle(objD);
                obj = gridHandle(obj);
                obj = obj.currentForm2Base(rad2deg([stepX,stepY]));
                % Test here for full 4pi grid, and if not, add the missing
                % axis/pole
                if obj1.isGrid4pi && ~obj.isGrid4pi
                    % Initially test for for the azimuth angles in the angle
                    % specified grids
                    if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
                        % We are expecting a symmetric x-grid from the TrueView
                        % transformation earlier, so only check for {-180,180}
                        xmin = min(obj.x);
                        xmax = max(obj.x);
                        tol = mean(diff(unique(obj.x)));
                        if abs(xmin + pi) > tol/10
                            % Insert a -pi cut from pi
                            obj = copyAndInsertXcut(obj,pi,-pi,tol);
                        end
                        if abs(xmax - pi) > tol/10
                            % Insert a pi cut from -pi
                            obj = copyAndInsertXcut(obj,-pi,pi,tol);
                        end
                    else
                        warning(['A full 4pi grid is expected, but not achieved.  Fix not yet implemented for gridType: ',obj.gridType]);
                    end
                end
            end
        end
        
        function obj = getBOR1pattern(obj1)
            % Function that expands the input FarField pattern into its BOR
            % components, and returns a FarField object which only contains
            % the BOR1 components.  The output field has the same th
            % samples as the input field, but only the principle ph cuts
            
            tol = 10^(-obj1.nSigDig+1);
            assert(strcmp(obj1.gridType,'PhTh'),'getBOR1pattern only operates on PhTh grid patterns');
            assert(abs(max(obj1.x) - min(obj1.x)) - 2*pi < tol,'The ph cuts must span 2*pi for BOR1 expansion');
            assert(obj1.isGridUniform,'A plaid, monotonic, uniform grid is expected for BOR1 expansion');
            assert(strcmp(obj1.coorType,'spherical'),'getBOR1pattern only operates on spherical coorType');
            
            Nph = obj1.Nx;
            Nth = obj1.Ny;
            th_vect = unique(obj1.y);
            ph_vect = unique(obj1.x);
            
            % Calculate the DFT in ph to get the BOR components
            % STore th variation in columns and BOR components row-wise
            [An,Bn,Cn,Dn] = deal(zeros(floor(Nph - 1)/2+1,Nth));
            [A1th,B1th,C1th,D1th] = deal(zeros(Nth,obj1.Nf));
            [BOR1power] = deal(zeros(1,obj1.Nf));
            for ff = 1:obj1.Nf
%                 for nn = 0:floor((Nph - 1)/2)
                for nn = 0:1 % Just get what is required for speed - can slot the rest in if more modes are needed later
                    sin_vect = sin(nn*ph_vect);
                    cos_vect = cos(nn*ph_vect);
                    for tt = 1:Nth
                        Gth_vect = obj1.E1((0:(Nph-1))*Nth+tt,ff);
                        Gph_vect = obj1.E2((0:(Nph-1))*Nth+tt,ff);
                        An(nn+1,tt) = 2/Nph.*sum(Gth_vect(:).*sin_vect(:));
                        Bn(nn+1,tt) = 2/Nph.*sum(Gth_vect(:).*cos_vect(:));
                        Cn(nn+1,tt) = 2/Nph.*sum(Gph_vect(:).*cos_vect(:));
                        Dn(nn+1,tt) = -2/Nph.*sum(Gph_vect(:).*sin_vect(:));
                    end
                end
                % Take the second rows as the BOR1 component - the first
                % row is BOR0
                A1th(:,ff) = An(2,:).';
                B1th(:,ff) = Bn(2,:).';
                C1th(:,ff) = Cn(2,:).';
                D1th(:,ff) = Dn(2,:).';
                
                BOR1power_integrand = 1./(2.*obj1.eta0).*(abs(A1th(:,ff)).^2 + abs(B1th(:,ff)).^2 + abs(C1th(:,ff)).^2 + abs(D1th(:,ff)).^2).*sin(th_vect);
                BOR1power(ff) = pi.*integral1D(th_vect,BOR1power_integrand);
            end
            % Build a suitable FarField class
            % For y-pol: A1 -> Gth and D1 -> Gph
            % For x-pol: B1 -> Gth and C1 -> Gph
            [PH,TH] = meshgrid([0,pi/2],th_vect);
            Eth = [B1th;A1th];  % First element corresponds to ph = 0, and second to ph = pi/2
            Eph = [C1th;D1th];
            obj = FarField(PH(:),TH(:),Eth,Eph,0.*Eth,obj1.freq,BOR1power,obj1.radEff,'spherical',obj1.polType,'PhTh',obj1.freqUnit,obj1.slant);
            obj.symmetryBOR1 = true;
        end
        
        function obj = expandBOR1pattern(obj1,phStepDeg)
            % Expands a BOR1 pattern, typically generated by
            % FarField.getBOR1pattern, into a 2*pi ph span
            
            if nargin < 2, phStepDeg = 5; end
            
            assert(obj1.symmetryBOR1,'Input object not BOR1 symmetric')
            assert(strcmp(obj1.gridType,'PhTh'),'BOR1 patterns must be specified on a PhTh grid')
            assert(strcmp(obj1.coorType,'spherical'),'BOR1 patterns must be specified in a Ludwig3 coordinate system')
            assert(isequal(unique(obj1.x),[0;pi/2]),'Expect ph cuts only at 0 and pi/2')
            assert(obj1.isGridUniform,'A plaid, monotonic, uniform grid is expected for BOR1 field expansion');
            
            phStep = deg2rad(phStepDeg);
            Nph = 2*pi/phStep + 1;
            Nth = obj1.Ny;
            ph_vect = linspace(0,2*pi,Nph);
            th_vect = obj1.y(1:Nth);
            [PH,TH] = meshgrid(ph_vect,th_vect);
            [A1,B1,C1,D1] = obj1.getBOR1comps;
            [Eth,Eph] = deal(zeros(Nth*Nph,obj1.Nf));
            for ff = 1:obj1.Nf
                Gth = bsxfun(@times,sin(PH),A1(:,ff)) + bsxfun(@times,cos(PH),B1(:,ff));
                Gph = bsxfun(@times,cos(PH),C1(:,ff)) - bsxfun(@times,sin(PH),D1(:,ff));
                Eth(:,ff) = Gth(:);
                Eph(:,ff) = Gph(:);
            end
            obj = FarField(PH(:),TH(:),Eth,Eph,0.*Eth,obj1.freq,obj1.Prad,obj1.radEff,'spherical',obj1.polType,'PhTh',obj1.freqUnit,obj1.slant);
        end
        
        function [A1,B1,C1,D1] = getBOR1comps(obj1)
            assert(obj1.symmetryBOR1,'Input object not BOR1 symmetric')
            assert(strcmp(obj1.gridType,'PhTh'),'BOR1 patterns must be specified on a PhTh grid')
            assert(strcmp(obj1.coorType,'spherical'),'BOR1 patterns must be specified in a Ludwig3 coordinate system')
            assert(isequal(unique(obj1.x),[0;pi/2]),'Expect ph cuts only at 0 and pi/2')
            assert(obj1.isGridUniform,'A plaid, monotonic, uniform grid is expected for BOR1 field expansion');
            Nth = obj1.Ny;
            A1 = obj1.E1(Nth+1:end,:);
            B1 = obj1.E1(1:Nth,:);
            C1 = obj1.E2(1:Nth,:);
            D1 = obj1.E2(Nth+1:end,:);
        end
        
        %% Format and other testers
        function y = isGridEqual(obj1,obj2)
            % Dont go for formal equality - floating point error just too much...
            tol = 10^(-obj1.nSigDig+1);
            if all(size(obj1.x) == size(obj2.x)) && all(size(obj1.y) == size(obj2.y))
                xEqual = abs(obj1.x - obj2.x) < tol;
                yEqual = abs(obj1.y - obj2.y) < tol;
                gridEqual = strcmp(obj1.gridType,obj2.gridType);
                fEqual = isequal(obj1.freq,obj2.freq);
                y = all(xEqual) && all(yEqual) && gridEqual && fEqual;
            else
                y = 0;
            end
        end
        
        function y = typesAreEqual(obj1,obj2)
            gridEqual = strcmp(obj1.gridType,obj2.gridType);
            coorEqual = strcmp(obj1.coorType,obj2.coorType);
            polEqual = strcmp(obj1.polType,obj2.polType);
            y = gridEqual && coorEqual && polEqual;
        end
        
        function y = isGrid4pi(obj)
            % Set to the PhTh coordinate system - this is how most data
            % will be generated anyway.
            % Very quick check - necessary but not always sufficient
            phRange = max(obj.phBase) - min(obj.phBase);
            thRange = max(obj.thBase) - min(obj.thBase);
            eps = 1e-4;
            y = ((abs(round(rad2deg(phRange)) - (360/2^(sum(abs([obj.symXZ,obj.symYZ]))))) < eps) & (abs(round(rad2deg(thRange)) - 180/2^abs(obj.symXY)) < eps)) |...
                ((abs(round(rad2deg(phRange)) - 180) < eps) & (abs(round(rad2deg(thRange)) - 360) < eps));
        end
        
        function y = isGridUniform(obj)
            % Test for a plaid, monotonic, uniform grid format
            
            if obj.Nx*obj.Ny == obj.Nang
                X = reshape(obj.x,obj.Ny,obj.Nx);
                Y = reshape(obj.y,obj.Ny,obj.Nx);
                % Test for equal rows in X, and equal columns in Y (plaid)
                rowEq = all(all(bsxfun(@eq,X,X(1,:))));
                colEq = all(all(bsxfun(@eq,Y,Y(:,1))));
                % Test for monotonic
                diffX = diff(X(1,:));
                diffY = diff(Y(:,1));
                monX = all(diffX>0);
                monY = all(diffY>0);
                % Test for uniform
                tol = 10^(-obj.nSigDig+1);
                unX = all(abs(diffX - median(diffX)) < tol);
                unY = all(abs(diffY - median(diffY)) < tol);
                % And all of them
                y = rowEq && colEq && monX && monY && unX && unY;
            else
                y = false;
            end
        end
        
        %% File Output methods
        writeGRASPcut(obj,pathName)
        
    end
    
    methods (Static = true)
        %% Farfield reading methods
        obj = readGRASPgrd(pathName);
        obj = readFEKOffe(pathName);
        obj = readCSTffs(pathName);
        obj = readGRASPcut(pathName,nr_freq,nr_cuts);
        obj = farFieldFromPowerPattern(x,y,P,freq,Pdim,coorType,polType,gridType,freqUnit,slant);
        
    end
    
    %% Internal helper functions
    methods (Access = private)
        
        function obj = setBase(obj)
            obj.xBase = obj.x;
            obj.yBase = obj.y;
            obj.NxBase = obj.Nx;
            obj.NyBase = obj.Ny;
            obj.phBase = obj.ph;
            obj.thBase = obj.th;
            obj.E1Base = obj.E1;
            obj.E2Base = obj.E2;
            obj.E3Base = obj.E3;
            obj.gridTypeBase = obj.gridType;
            obj.coorTypeBase = obj.coorType;
            obj.polTypeBase = obj.polType;
        end
        
        function [ph,th] = getPhThCurrent(obj)
            % This does not operate on the base grid, but instead the
            % current grid
            if strcmp(obj.gridType,'PhTh')
                ph = obj.x;
                th = obj.y;
            else
                handle2DirCos = str2func([obj.gridType,'2DirCos']);
                [u,v,w] = handle2DirCos(obj.x,obj.y);
                [ph,th] = DirCos2PhTh(u,v,w);
            end
        end
        
        % Set the names of the 2 grid components
        function [xname,yname] = setXYnames(obj)
            switch obj.gridType
                case 'PhTh'
                    xname = '\phi';
                    yname = '\theta';
                case 'DirCos'
                    xname = 'u';
                    yname = 'v';
                case 'AzEl'
                    xname = 'az';
                    yname = 'el';
                case 'ElAz'
                    xname = '\epsilon';
                    yname = '\alpha';
                case 'TrueView'
                    xname = 'Xg=\theta cos(\phi)';
                    yname = 'Yg=\theta sin(\phi)';
                case 'ArcSin'
                    xname = 'Xg=asin(u)';
                    yname = 'Yg=asin(v)';
            end
        end
        
        % Set the names of the 2 farfield components based on the
        % polarization type.  Names used for info and plotting.
        function [E1name,E2name] = setEnames(obj)
            switch obj.polType
                case 'circular'
                    E1name = 'Elh';
                    E2name = 'Erh';
                case 'slant'
                    E1name = 'Exp';
                    E2name = 'Eco';
                case 'linear'
                    switch obj.coorType
                        case 'spherical'
                            E1name = 'Eth';
                            E2name = 'Eph';
                        case 'Ludwig1'
                            E1name = 'Ex';
                            E2name = 'Ey';
                        case 'Ludwig2AE'
                            E1name = 'Eaz';
                            E2name = 'Eel';
                        case 'Ludwig2EA'
                            E1name = 'Eal';
                            E2name = 'Eep';
                        case 'Ludwig3'
                            E1name = 'Eh';
                            E2name = 'Ev';
                    end
                    
                otherwise
                    error(['Unknown polType property: ', obj.polType]);
            end
        end
        
        function [xRangeType,yRangeType] = setRangeTypes(obj)
            % Try to figure out what the current rangeType is.
            % Not much error checking is done - assume somewhat
            % sensible inputs are provided most of the time.
            xRangeType = 'sym';
            if (strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz'))
                if min(obj.x) >= 0
                    xRangeType = 'pos';
                end
                if max(obj.y) - min(obj.y) <= pi+median(diff(unique(obj.y)))/2
                    yRangeType = '180';
                else
                    yRangeType = '360';
                end
            else
                yRangeType = [];
            end
        end
        
    end
    
end