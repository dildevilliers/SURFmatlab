classdef FarField
    % Based largely on information in CST help file: Farfield Calculation
    % Overview, and the AMTA paper in this folder 'COORDINATE SYSTEM
    % PLOTTING FOR ANTENNA MEASUREMENTS', GF Masters and SF Gregson
    % Currently assumes direction of propagation is the global z-axis for 
    % all directional polarization types.
    
    properties
        x(:,1) double {mustBeReal, mustBeFinite}
        y(:,1) double {mustBeReal, mustBeFinite}
        r(1,1) double {mustBeReal, mustBeFinite} = 1
        E1(:,:) double %{mustBeFinite}
        E2(:,:) double %{mustBeFinite}
        E3(:,:) double %{mustBeFinite}
        freq(1,:) double {mustBeReal, mustBeFinite} = 1
        Prad(1,:) double {mustBeReal, mustBeFinite} = 4*pi
        radEff(1,:) double {mustBeReal, mustBeFinite} = 1
        coorSys(1,:) char {mustBeMember(coorSys,{'spherical','Ludwig1','Ludwig2AE','Ludwig2EA','Ludwig3'})} = 'spherical'
        polType(1,:) char {mustBeMember(polType,{'linear','circular','slant'})} = 'linear'
        slant(1,1) double {mustBeReal, mustBeFinite} = 0   % slant angle in radians - measured between Exp and E1
        gridType(1,:) char {mustBeMember(gridType,{'PhTh','DirCos','AzEl','ElAz','TrueView','ArcSin'})} = 'PhTh'    
        freqUnit(1,:) char {mustBeMember(freqUnit,{'Hz','kHz','MHz','GHz','THz'})} = 'Hz'    
    end
    
    properties (SetAccess = private)
        E1name  % ['Eth', 'Ex', 'Eaz', 'Eal', 'Eh', 'Elh', 'Exp'] - depends on coorSys and polType
        E2name  % ['Eph', 'Ey', 'Eel', 'Eep', 'Ev', 'Erh', 'Eco'] - depends on coorSys and polType
        xname
        yname
        ph
        th
        freqHz
        Nf      % number of frequencies
        Nx     % number of unique x angles
        Ny     % number of unique y angles
        Nang    % number of angle combinations [Nx x Ny]
        radEff_dB       % radiation efficiency in dB [1 x Nf]
        Directivity_dBi % directivity in dBi [1 x Nf]
        Gain_dB         % Gain in dB [1 x Nf]
        xRangeType     % 'sym' or 'pos' 
        yRangeType     % '180' or '360' 

    end
    
    properties (SetAccess = private, Hidden = true)
        % Keep the input data here to not lose some info when going through
        % a DirCos projection and back...
        xBase   
        yBase 
        gridTypeBase
        E1Base
        E2Base
        E3Base
        coorSysBase
        polTypeBase
    end
        
        properties (Constant = true, Hidden = true)
        c0 = physconst('Lightspeed');
        eps0 = 8.854187817000001e-12;
        mu0 = 1.256637061435917e-06;
        eta0 = 3.767303134749689e+02;
    end
    
    methods
        % Make a basic constructor method 
        function obj = FarField(x,y,E1,E2,E3,freq,Prad,radEff)
            
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
            else
                Nang = Nang_x;
            end
            if ~isequal(Nf_E1,Nf_E2,Nf_E3)
                error('E1, E2, and E3 must have the same number of columns');
            end
            Nf = Nf_E1;
            if nargin < 6
                freq = ones(1,Nf).*obj.freq;
                warning('freq not specified - using default for all columns');
            else
                if Nf ~= length(freq(1,:))
                    error('freq must have the same number of columns as E1 and E2');
                end
            end
            if nargin < 7
                Prad = ones(1,Nf).*obj.Prad;
                warning('Prad not specified - using default for all columns');
            else
                if Nf ~= length(Prad(1,:))
                    error('Prad must have the same number of columns as E1 and E2');
                end
            end
            if nargin < 8
                radEff = ones(1,Nf).*obj.radEff;
                warning('radEff not specified - using default for all columns');
            else
                if Nf ~= length(radEff(1,:))
                    error('radEff must have the same number of columns as E1 and E2');
                end
            end
            if Nf_x > 1 || Nf_y > 1
                warning('Only using first column of th and ph since they must be equal for all frequencies anyway');
                x = x(:,1);
                y = y(:,1);
            end
            obj.x = x;
            obj.y = y;
%             [obj.th, obj.ph] = obj.getphth();
            obj.E1 = E1;
            obj.E2 = E2;
            obj.E3 = E3;
            obj.freq = freq;
            obj.Prad = Prad;
            obj.radEff = radEff;
            obj.radEff_dB = dB10(radEff);
            obj.Nf = Nf;
            obj.Nang = Nang;
            obj.Nx = length(unique(x));
            obj.Ny = length(unique(y));
            obj.Directivity_dBi = dB10(max(obj.getDirectivity()));
            obj.Gain_dB = dB10(max(obj.getGain()));
            obj = setEnames(obj);
            obj = setXYnames(obj);
            obj = setBase(obj);
            obj = setFreq(obj);
            obj = setPhTh(obj);
            obj = setRangeTypes(obj);
        end
        
        
        %% Pattern getters
        function [E1field, E2field, E3field] = getEfield(obj)
            % function [E1field, E2field, E3field] = getEfield(obj)
            % Returns the Efield matrices of size [Nang x Nf]
            % Efield = E*exp(-jkr)/r
            k = 2.*pi.*obj.freq./obj.c0;
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
            handle2DirCos = str2func(['FarField.',obj.gridTypeBase,'2DirCos']);
            [u,v,w] = handle2DirCos(obj.xBase,obj.yBase);
        end
        
        function [ph, th] = getPhTh(obj)
            switch obj.gridTypeBase
                case 'PhTh'
                    ph = obj.xBase;
                    th = obj.yBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [ph,th] = FarField.DirCos2PhTh(u,v,w);
            end
        end
        
        function [az, el] = getAzEl(obj)
            switch obj.gridTypeBase
                case 'AzEl'
                    el = obj.yBase;
                    az = obj.xBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [az,el] = FarField.DirCos2AzEl(u,v,w);
            end
        end
        
        function [ep, al] = getElAz(obj)
            switch obj.gridTypeBase
                case 'ElAz'
                    ep = obj.xBase;
                    al = obj.yBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [ep,al] = FarField.DirCos2ElAz(u,v,w);
            end
        end
        
        function [Xg, Yg] = getTrueView(obj)
            switch obj.gridTypeBase
                case 'TrueView'
                    Xg = obj.xBase;
                    Yg = obj.xBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [Xg,Yg] = FarField.DirCos2TrueView(u,v,w);
            end
        end
        
        function [asinu, asinv] = getArcSin(obj)
            switch obj.gridTypeBase
                case 'ArcSin'
                    asinu = obj.xBase;
                    asinv = obj.yBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [asinu,asinv] = FarField.DirCos2ArcSin(u,v,w);
            end
        end
        
        %% Grid transformation setters
        function obj = grid2PhTh(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'PhTh')
                [obj.x,obj.y] = getPhTh(obj);
                obj.gridType = 'PhTh';
                obj = setXYnames(obj);
                obj = setRangeTypes(obj);
            end
        end
        
        function obj = grid2DirCos(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'DirCos')
                [obj.x,obj.y] = getDirCos(obj);
                obj.gridType = 'DirCos';
                obj = setXYnames(obj);
                obj = setRangeTypes(obj);
            end
        end
        
        function obj = grid2AzEl(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'AzEl')
                [obj.x,obj.y] = getAzEl(obj);
                obj.gridType = 'AzEl';
                obj = setXYnames(obj);
                obj = setRangeTypes(obj);
            end
        end
        
        function obj = grid2ElAz(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'ElAz')
                [obj.x,obj.y] = getElAz(obj);
                obj.gridType = 'ElAz';
                obj = setXYnames(obj);
                obj = setRangeTypes(obj);
            end
        end
        
        function obj = grid2TrueView(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'TrueView')
                [obj.x,obj.y] = getTrueView(obj);
                obj.gridType = 'TrueView';
                obj = setXYnames(obj);
                obj = setRangeTypes(obj);
            end
        end
        
        function obj = grid2ArcSin(obj)
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'ArcSin')
                [obj.x,obj.y] = getArcSin(obj);
                obj.gridType = 'ArcSin';
                obj = setXYnames(obj);
                obj = setRangeTypes(obj);
            end
        end
        
        %% Grid range shifters
        function obj = sortGrid(obj)
            obj = roundGrid(obj);
            [~,iSort] = sortrows([obj.x,obj.y],[1 2]);
            obj.x = obj.x(iSort);
            obj.y = obj.y(iSort);
            obj.Nx = numel(unique(obj.x));
            obj.Ny = numel(unique(obj.y));
            obj.E1 = obj.E1(iSort);
            obj.E2 = obj.E2(iSort);
            obj.E3 = obj.E3(iSort);
        end
        
        function obj = roundGrid(obj,nSigDig)
            % Round to some significant digits for sorting (some issues can
            % arise in deg2rad and rad2deg
            if nargin < 2
                nSigDig = 10;
            end
            xRound = round(obj.x*10^nSigDig)/10^nSigDig;
            yRound = round(obj.y*10^nSigDig)/10^nSigDig;
            obj.x = xRound;
            obj.y = yRound;
        end
        
        function obj = setRangeTypes(obj)
            % Try to figure out what the current rangeType is.
            % Not much error checking is done - assume somewhat
            % sensible inputs are provided most of the time.
            obj.xRangeType = 'sym';
            if (strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz'))
                if min(obj.x) >= 0
                    obj.xRangeType = 'pos';
                end
                if max(obj.y) - min(obj.y) <= pi
                    obj.yRangeType = '180';
                else
                    obj.yRangeType = '360';
                end
            else
                obj.yRangeType = [];
            end
        end
        
        obj = setXrange(obj,type)
        obj = setYrange(obj,type)
        
        %% Coordinate system getters
        function [Eth, Eph, Er] = getEspherical(obj)
            [Ph,Th] = getPhTh(obj);
            TH = repmat(Th(:,1),1,obj.Nf);
            PH = repmat(Ph(:,1),1,obj.Nf);
            % Change to the Base values here....
            switch obj.coorSysBase
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
            switch obj.coorSysBase
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
            switch obj.coorSysBase
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
                    % Sort out singularities poles
                    phPoles = deg2rad([-270,-90,90,270].');
                    poleMat = [ones(4,1).*deg2rad(90),phPoles]; % [th=90,ph]
                    [~,iPole] = ismember(poleMat,[obj.th,obj.ph],'rows');
                    iPole = iPole(iPole>0);
                    [Eaz(iPole,:),Eel(iPole,:),E3(iPole,:)] = deal(0);
%                     Eaz(iPole,:) = Eth(iPole,:) + Eph(iPole,:);
            end
        end
        
        function [Eal, Eep, E3] = getELudwig2EA(obj)
            switch obj.coorSysBase
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
                    % Sort out singularities poles
                    phPoles = deg2rad([-360,-180,0,180,360].');
                    poleMat = [ones(5,1).*deg2rad(90),phPoles]; % [th=90,ph]
                    [~,iPole] = ismember(poleMat,[obj.th,obj.ph],'rows');
                    iPole = iPole(iPole>0);
                    [Eal(iPole,:),Eep(iPole,:),E3(iPole,:)] = deal(0);
            end
        end
        
        function [Eh, Ev, E3] = getELudwig3(obj)
            switch obj.coorSysBase
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
            if ~strcmp(obj.coorSys,'spherical')
                [obj.E1,obj.E2,obj.E3] = getEspherical(obj);
                obj.E3 = zeros(size(obj.E1));
                obj.coorSys = 'spherical';
                obj = setEnames(obj);
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        function obj = coor2Ludwig1(obj,setStdGrid)
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorSys,'Ludwig1')
                [obj.E1,obj.E2,obj.E3] = getELudwig1(obj);
                obj.coorSys = 'Ludwig1';
                obj = setEnames(obj);
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        function obj = coor2Ludwig2AE(obj,setStdGrid)
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorSys,'Ludwig2AE')
                [obj.E1,obj.E2,obj.E3] = getELudwig2AE(obj);
                obj.coorSys = 'Ludwig2AE';
                obj = setEnames(obj);
            end
            if setStdGrid
                obj = obj.grid2AzEl;
            end
        end
        
        function obj = coor2Ludwig2EA(obj,setStdGrid)
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorSys,'Ludwig2EA')
                [obj.E1,obj.E2,obj.E3] = getELudwig2EA(obj);
                obj.coorSys = 'Ludwig2EA';
                obj = setEnames(obj);
            end
            if setStdGrid
                obj = obj.grid2ElAz;
            end
        end
        
        function obj = coor2Ludwig3(obj,setStdGrid)
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorSys,'Ludwig3')
                [obj.E1,obj.E2,obj.E3] = getELudwig3(obj);
                obj.coorSys = 'Ludwig3';
                obj = setEnames(obj);
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        %% Polarization type getters
        function [E1lin, E2lin, E3lin] = getElin(obj)
            switch obj.polTypeBase
                case 'linear'
                    E1lin = obj.E1Base;
                    E2lin = obj.E2Base;
                case 'circular'
                    Del = 2*1i;
                    E1lin = sqrt(2)./Del.*(1i.*obj.E1Base + 1i.*obj.E2Base);
                    E2lin = sqrt(2)./Del.*(-obj.E1Base + obj.E2Base);
                case 'slant'
                    PSI = ones(size(obj.E1Base)).*obj.slant;
                    Del = 1;
                    E1lin = 1./Del.*(cos(PSI).*obj.E1Base + sin(PSI).*obj.E2Base);
                    E2lin = 1./Del.*(-sin(PSI).*obj.E1Base + cos(PSI).*obj.E2Base);
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
                obj = setEnames(obj);
            end
        end
        
        function obj = pol2circular(obj)
            if ~strcmp(obj.polType,'circular')
                [obj.E1,obj.E2,obj.E3] = getEcircular(obj);
                obj.polType = 'circular';
                obj = setEnames(obj);
            end
        end
        
        function obj = pol2slant(obj)
            if ~strcmp(obj.polType,'slant')
                [obj.E1,obj.E2,obj.E3] = getEslant(obj);
                obj.polType = 'slant';
                obj = setEnames(obj);
            end
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
            obj.coorSys = obj.coorSysBase;
            obj.polType = obj.polTypeBase;
            obj = setEnames(obj);
            obj = setXYnames(obj);
            obj = setRangeTypes(obj);
        end
        
        function obj = grid2Base(obj)
            % Evaluate the current object (pol and coor) on the base grid
            coorSysIn = obj.coorSys;
            polTypeIn = obj.polType;
            coorSysH = str2func(['coor2',coorSysIn]);
            polTypeH = str2func(['pol2',polTypeIn]);
            obj = obj.reset2Base;
            % Keep the current coorSys and polType
            obj = coorSysH(obj,false);
            obj = polTypeH(obj);
            obj = setEnames(obj);
        end
        
        function objNew = currentForm2Base(obj,stepDeg,hemisphere)
            % Sets the base to the current format. Resamples the field on a
            % regular plaid grid, and makes this the new base grid. This is
            % typically then not where actual samples where, but instead
            % interpolated values.
            nSigDig = 10;
            if nargin < 2
                % Get a step from the current object
                if strcmp(obj.gridTypeBase,'DirCos') || strcmp(obj.gridType,'ArcSin')
                    stepX = asin(min(abs(diff(unique(obj.xBase)))));
                    stepY = asin(min(abs(diff(unique(obj.yBase)))));
                else
                    % Sort out rounding errors for degrees
                    stepX = deg2rad(round(rad2deg(min(abs(diff(unique(obj.xBase)))))*10^nSigDig)/10^nSigDig);
                    stepY = deg2rad(round(rad2deg(min(abs(diff(unique(obj.yBase)))))*10^nSigDig)/10^nSigDig);
                end
                hemisphere = 'top';
            elseif nargin < 3
                hemisphere = 'top';
                if numel(stepDeg) == 1
                    [stepX,stepY] = deal(deg2rad(stepDeg));
                elseif numel(stepDeg) == 2
                    stepX = deg2rad(stepDeg(1));
                    stepY = deg2rad(stepDeg(2));
                end
            end
            
            if strcmp(obj.gridType,'DirCos') || strcmp(obj.gridType,'ArcSin')
                stepX = sin(stepX);
                stepY = sin(stepY);
            end
            
            % Build the new grid
            xivect = min(obj.x):stepX:max(obj.x);
            yivect = min(obj.y):stepY:max(obj.y);
            [Xi,Yi] = meshgrid(xivect,yivect);
            xi = Xi(:);
            yi = Yi(:);
            Nxi = numel(xivect);
            Nyi = numel(yivect);
            % Interpolate the fields
            [E1grid,E2grid,E3grid] = deal(zeros(Nxi*Nyi,obj.Nf));
            for ff = 1:obj.Nf
                E1grid(:,ff) = interpolateGrid(obj,'E1',xi,yi,ff,hemisphere);
                E2grid(:,ff) = interpolateGrid(obj,'E2',xi,yi,ff,hemisphere);
                E3grid(:,ff) = interpolateGrid(obj,'E3',xi,yi,ff,hemisphere);
            end
            % Populate the new farField object
            objNew = obj;
            objNew.x = xi(:);
            objNew.y = yi(:);
            objNew.Nx = Nxi;
            objNew.Ny = Nyi;
            objNew.Nang = Nxi*Nyi;
            objNew.E1 = E1grid;
            objNew.E2 = E2grid;
            objNew.E3 = E3grid;
            objNew = setBase(objNew);
            % Update the current form to the base form
            objNew = reset2Base(objNew);
        end
        
        
        %% Plotting methods
        plot(obj,varargin)
        plotJones(obj1,obj2,varargin)  % Much to do here still...
        
        %% Interpolation methods
%         [Z] = interpolateGrid(obj,xi,yi,gridType,output,freqIndex)
        [Z] = interpolateGrid(obj,output,xi,yi,varargin)
        
        %% Maths
        function obj = plus(obj1,obj2)
            obj1 = resetToBase(obj1);
            obj2 = resetToBase(obj2);
            
            if isGridEqual(obj1,obj2)
                obj = obj1;
                obj.E1 = obj1.E1 + obj2.E1;
                obj.E2 = obj1.E2 + obj2.E2;
                obj.E3 = obj1.E3 + obj2.E3;
                obj.Prad = obj1.Prad + obj2.Prad;
                Pt = obj1.Prad./obj1.radEff + obj2.Prad./obj2.radEff;
                obj.radEff = obj.Prad./Pt;
                obj.radEff_dB = dB10(obj.radEff);
                obj.Directivity_dBi = dB10(max(obj.getDirectivity()));
                obj.Gain_dB = dB10(max(obj.getGain()));
            else 
                error('Can only add FarFields with equal base grids')
            end
        end
        
        function obj = times(obj1,obj2)
            obj1 = resetToBase(obj1);
            obj2 = resetToBase(obj2);
            
            if isGridEqual(obj1,obj2)
                obj = obj1;
                obj.E1 = obj1.E1.*obj2.E1;
                obj.E2 = obj1.E2.*obj2.E2;
                obj.E3 = obj1.E3.*obj2.E3;
                % ToDo - calculate power through field integration
                % For now just normalise to 4pi
                obj.Prad = ones(size(obj.Prad)).*4*pi;
                obj.radEff = ones(size(obj.Prad));
                obj.radEff_dB = dB10(obj.radEff);
                obj.Directivity_dBi = dB10(max(obj.getDirectivity()));
                obj.Gain_dB = dB10(max(obj.getGain()));
            else 
                error('Can only multiply FarFields with equal base grids')
            end
        end
        
        %% Format and other testers
        function y = isGridEqual(obj1,obj2)
            xEqual = isequal(obj1.x,obj2.x);
            yEqual = isequal(obj1.y,obj2.y);
            gridEqual = strcmp(obj1.gridType,obj2.gridType);
            fEqual = isequal(obj1.freq,obj2.freq);
            y = xEqual && yEqual && gridEqual && fEqual;
        end
        
    end
    
    methods (Static = true)
       %% Farfield reading methods
        obj = readGRASPgrd(pathName); 
        
        %% Coordinate transformers
        function [u,v,w] = PhTh2DirCos(ph,th)
            u = real(sin(th).*cos(ph));
            v = real(sin(th).*sin(ph));
            w = real(cos(th));
        end
        
        function [u,v,w] = DirCos2DirCos(u,v,w)
            if nargin < 3
                w = sqrt(1 - u.^2 - v.^2);
            end
        end
        
        function [u,v,w] = AzEl2DirCos(az,el)
            u = real(sin(az).*cos(el));
            v = real(sin(el));
            w = real(cos(az).*cos(el));
        end
        
        function [u,v,w] = ElAz2DirCos(ep,al)
            u = real(sin(al));
            v = real(cos(al).*sin(ep));
            w = real(cos(al).*cos(ep));
        end
        
        function [u,v,w] = TrueView2DirCos(x,y)
            Th = sqrt(x.^2 + y.^2);
            Ph = atan2(y,x);
            u = real(sin(Th).*cos(Ph));
            v = real(sin(Th).*sin(Ph));
            w = real(cos(Th));
        end
        
        function [u,v,w] = ArcSin2DirCos(x,y)
            u = sin(x);
            v = sin(y);
            w = real(sqrt(1 - (sin(x).^2 + sin(y).^2)));
        end
        
        function [ph,th] = DirCos2PhTh(u,v,w)
            ph = real(atan2(v,u));
            th = real(acos(w));
        end
        
        function [az,el] = DirCos2AzEl(u,v,w)
            el = real(asin(v));
            az = real(atan2(u,w));
        end
        
        function [ep,al] = DirCos2ElAz(u,v,w)
            al = real(asin(u));
            ep = real(atan2(v,w));
        end
        
        function [Xg,Yg] = DirCos2TrueView(u,v,w)
            Ph = atan2(v,u);
            Th = acos(w);
            Xg = real(Th.*cos(Ph));
            Yg = real(Th.*sin(Ph));
        end

        function [asinu,asinv] = DirCos2ArcSin(u,v,w)
            asinu = real(asin(u));
            asinv = real(asin(v));
        end
    end
    
    %% Internal helper functions
    methods (Access = private)
        
        function obj = setBase(obj)
            obj.xBase = obj.x;
            obj.yBase = obj.y;
            obj.gridTypeBase = obj.gridType;
            obj.E1Base = obj.E1;
            obj.E2Base = obj.E2;
            obj.E3Base = obj.E3;
            obj.coorSysBase = obj.coorSys;
            obj.polTypeBase = obj.polType;
        end
        
        function obj = setFreq(obj)
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
            obj.freqHz = obj.freq*freqMult;
        end
        
        function obj = setPhTh(obj)
            [obj.ph, obj.th] = getPhTh(obj);
        end
        
        % Set the names of the 2 grid components 
        function obj = setXYnames(obj)
            switch obj.gridType
                case 'PhTh'
                    obj.xname = '\phi';
                    obj.yname = '\theta';
                case 'DirCos'
                    obj.xname = 'u';
                    obj.yname = 'v';
                case 'AzEl'
                    obj.xname = 'az';
                    obj.yname = 'el';
                case 'ElAz'
                    obj.xname = '\epsilon';
                    obj.yname = '\alpha';
                case 'TrueView'
                    obj.xname = 'Xg=\theta cos(\phi)';
                    obj.yname = 'Yg=\theta sin(\phi)';
                case 'ArcSin'
                    obj.xname = 'Xg=asin(u)';
                    obj.yname = 'Yg=asin(v)';
            end
        end
        
        % Set the names of the 2 farfield components based on the
        % polarization type.  Names used for info and plotting.
        function obj = setEnames(obj)
            switch obj.polType
                case 'circular'
                    obj.E1name = 'Elh';
                    obj.E2name = 'Erh';
                case 'slant'
                    obj.E1name = 'Exp';
                    obj.E2name = 'Eco';
                case 'linear'
                    switch obj.coorSys
                        case 'spherical'
                            obj.E1name = 'Eth';
                            obj.E2name = 'Eph';
                        case 'Ludwig1'
                            obj.E1name = 'Ex';
                            obj.E2name = 'Ey';
                        case 'Ludwig2AE'
                            obj.E1name = 'Eaz';
                            obj.E2name = 'Eel';
                        case 'Ludwig2EA'
                            obj.E1name = 'Eal';
                            obj.E2name = 'Eep';
                        case 'Ludwig3'
                            obj.E1name = 'Eh';
                            obj.E2name = 'Ev';
                    end
                    
                otherwise
                    error(['Unknown polType property: ', obj.polType]);
            end
        end
        
    end
    
end