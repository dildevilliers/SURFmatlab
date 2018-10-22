classdef FarField
    % Based largely on information in CST help file: Farfield Calculation
    % Overview, and the AMTA paper in this folder 'COORDINATE SYSTEM
    % PLOTTING FOR ANTENNA MEASUREMENTS', GF Masters and SF Gregson
    % Currently assumes direction of propagation is the global z-axis for 
    % all directional polarization types.
    
    properties
        th      % th angle vector [Nang x 1] - typically in [0,pi] rad
        ph      % ph angle vector [Nang x 1] - typically in [0,2*pi] rad
        r = 1;  % radius for E-field calculation in m  
        E1      % first component of E-field pattern [Nang x Nf] - depends on polBase and polType
        E2      % second component of E-field pattern [Nang x Nf] - depends on polBase and polType
        E3      % third component of E-field pattern [Nang x Nf] - only exists for polBase = Ludwig1
        freq = 1;               % frequency vector in Hz [1 x Nf]
        Prad = 4*pi;            % radiated power vector in W [1 x Nf]
        radEff = 1;             % antenna radiation efficiency pu [1 x Nf]
        polBase = 'spherical';  % polarization base of field [('spherical'), 'Ludwig1', 'Ludwig2AE', 'Ludwig2EA', 'Ludwig3']
        polType = 'linear';     % polarization type of the field [('linear'), 'circular', 'slant']
        slant   = 0;            % slant angle in radians - measured between Exp and E1
    end
    
    properties (SetAccess = private)
        E1name  % ['Eth', 'Ex', 'Eaz', 'Eal', 'Eh', 'Elh', 'Exp'] - depends on polBase and polType
        E2name  % ['Eph', 'Ey', 'Eel', 'Eep', 'Ev', 'Erh', 'Eco'] - depends on polBase and polType
%         E1field % E1 electric field strength - E1*exp(-jkr)/r
%         E2field % E2 electric field strength - E2*exp(-jkr)/r
%         E3field % E2 electric field strength - E2*exp(-jkr)/r
        Nf      % number of frequencies
        Nth     % number of unique theta angles
        Nph     % number of unique phi angles
        Nang    % number of angle combinations [Nth x Nph]
        radEff_dB       % radiation efficiency in dB [1 x Nf]
        Directivity_dBi % directivity in dBi [1 x Nf]
        Gain_dB         % Gain in dB [1 x Nf]
    end
    
    properties (Constant = true, Hidden = true)
        c0 = physconst('Lightspeed');
        eps0 = 8.854187817000001e-12;
        mu0 = 1.256637061435917e-06;
        eta0 = 3.767303134749689e+02;
    end
    
    methods
        % Make a basic constructor method 
        function obj = FarField(th,ph,E1,E2,E3,freq,Prad,radEff)
            
            % function obj = FarField(th,ph,E1,E2,E3,freq,Prad,radEff)
            % Constructor method for the FarField object
            % Required inputs:
            % th: column vector [Nang x 1] of th angles in rad
            % ph: column vector [Nang x 1] of ph angles in rad
            % E1: First E-field pattern component (exp(jkr)/r suppressed) of size [Nang x Nf]
            % E2: Second E-field pattern component (exp(jkr)/r suppressed) of size [Nang x Nf]
            %
            % Optional inputs - can be left empty or omitted
            % 'E3': Third E-field pattern component (exp(jkr)/r suppressed) of size [Nang x Nf]
            % 'freq': vector [1 x Nf] of frequencies where the fields are defined in Hz
            % 'Prad': vector [1 x Nf] of radiated powers in W
            % 'radEff': vector [1 x Nf] of radiation efficiencies
            
            % Basic input error checking
            [Nang_th, Nf_th] = size(th);
            [Nang_ph, Nf_ph] = size(ph);
            [Nang_E1, Nf_E1] = size(E1);
            [Nang_E2, Nf_E2] = size(E2);
            if nargin < 5 || isempty(E3)
                E3 = zeros(size(E2));
            end
            [Nang_E3, Nf_E3] = size(E3);
            if ~isequal(Nang_th,Nang_ph,Nang_E1,Nang_E2,Nang_E3)
                error('All inputs must have the same number of rows');
            else
                Nang = Nang_th;
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
            if Nf_th > 1 || Nf_ph > 1
                warning('Only using first column of th and ph since they must be equal for all frequencies anyway');
                th = th(:,1);
                ph = ph(:,1);
            end
            obj.th = th;
            obj.ph = ph;
            obj.E1 = E1;
            obj.E2 = E2;
            obj.E3 = E3;
            obj.freq = freq;
            obj.Prad = Prad;
            obj.radEff = radEff;
            obj.radEff_dB = dB10(radEff);
            obj.Nf = Nf;
            obj.Nang = Nang;
            obj.Nth = length(unique(th));
            obj.Nph = length(unique(ph));
            obj.Directivity_dBi = dB10(max(obj.getDirectivity()));
            obj.Gain_dB = dB10(max(obj.getGain()));
            obj = setEnames(obj);
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
        
        function [AR,ARinv] = getAxialRatio(obj)
            % function [AR, ARinv] = getAxialRatio(obj)
            % returns the Axial Ratio (linear) in AR and the inverted Axial Ratio in ARinv [Nang x Nf]
           AR = sqrt((abs(obj.E1).^2 + abs(obj.E2).^2 + abs(obj.E1.^2 + obj.E2.^2))./(abs(obj.E1).^2 + abs(obj.E2).^2 - abs(obj.E1.^2 + obj.E2.^2)));
           ARinv = sqrt((abs(obj.E1).^2 + abs(obj.E2).^2 - abs(obj.E1.^2 + obj.E2.^2))./(abs(obj.E1).^2 + abs(obj.E2).^2 + abs(obj.E1.^2 + obj.E2.^2)));
        end
        
        %% Angle and projection getters
        function [u, v, w] = getUVW(FF)
            u = sin(FF.th).*cos(FF.ph);
            v = sin(FF.th).*sin(FF.ph);
            w = cos(FF.th);
        end
        
        function [el, az] = getElAz(FF)
            [u,v,w] = getUVW(FF);
            el = asin(v);
            az = atan2(u,w);
        end
        
        function [ep, al] = getEpAl(FF)
            [u,v,w] = getUVW(FF);
            al = asin(u);
            ep = atan2(v,w);
        end
        
        function [Xg, Yg] = getXgYg(FF)
            Xg = FF.th.*cos(FF.ph);
            Yg = FF.th.*sin(FF.ph);
        end
        
        %% Polarization transformation methods
        % All take the current polarization base of the input object and 
        % transforms it to the required base
        % Strategy is to make the spherical components first, and then go
        % to whatever the requested version is...
        % First change of polBase
        function [Eth, Eph] = getEthEph(obj)
            TH = repmat(obj.th(:,1),1,obj.Nf);
            PH = repmat(obj.ph(:,1),1,obj.Nf);
            switch obj.polBase
                case 'spherical'
                    Eth = obj.E1;
                    Eph = obj.E2;
                case 'Ludwig1'
                    Eth = cos(TH).*cos(PH).*obj.E1 + cos(TH).*sin(PH).*obj.E2 - sin(TH).*obj.E3;
                    Eph = -sin(PH).*obj.E1 + cos(PH).*obj.E2;
                case 'Ludwig2AE'
                    cosEl = sqrt(1 - sin(TH).^2.*sin(PH).^2);
                    Del = cos(PH).^2 + cos(TH).^2.*sin(PH).^2;
                    Eth = (cosEl./Del).*(cos(PH).*obj.E1 + cos(TH).*sin(PH).*obj.E2);
                    Eph = (cosEl./Del).*(-cos(TH).*sin(PH).*obj.E1 + cos(PH).*obj.E2);
                case 'Ludwig2EA'
                    cosAl = sqrt(1 - sin(TH).^2.*cos(PH).^2);
                    Del = cos(TH).^2.*cos(PH).^2 + sin(PH).^2;
                    Eth = (cosAl./Del).*(cos(TH).*cos(PH).*obj.E1 + sin(PH).*obj.E2);
                    Eph = (cosAl./Del).*(-sin(PH).*obj.E1 + cos(TH).*cos(PH).*obj.E2);
                case 'Ludwig3'
                    Del = 1;
                    Eth = (1./Del).*(cos(PH).*obj.E1 + sin(PH).*obj.E2);
                    Eph = (1./Del).*(-sin(PH).*obj.E1 + cos(PH).*obj.E2);
                otherwise
                    error(['Unknown polBase property: ', obj.polBase]);
            end
        end
        
        function obj = pol2spherical(obj)
            if ~strcmp(obj.polType,'spherical')
                [Eth, Eph] = getEthEph(obj);
                obj.E1 = Eth;
                obj.E2 = Eph;
                obj.E3 = zeros(size(obj.E1));
                obj.polBase = 'spherical';
                obj = setEnames(obj);
            end
        end
        
        function obj = pol2Ludwig1(obj)
            if ~strcmp(obj.polType,'Ludwig1')
                [Eth, Eph] = getEthEph(obj);
                TH = repmat(obj.th(:,1),1,obj.Nf);
                PH = repmat(obj.ph(:,1),1,obj.Nf);
                % Assume farfield so no radial E-field
                obj.E1 = cos(TH).*cos(PH).*Eth - sin(PH).*Eph;
                obj.E2 = cos(TH).*sin(PH).*Eth + cos(PH).*Eph;
                obj.E3 = zeros(size(obj.E1));   % Strict definition in the Ludwig paper
                obj.polBase = 'Ludwig1';
                obj = setEnames(obj);
            end
        end
        
        function obj = pol2Ludwig2AE(obj)
            if ~strcmp(obj.polType,'Ludwig2AE')
                [Eth, Eph] = getEthEph(obj);
                TH = repmat(obj.th(:,1),1,obj.Nf);
                PH = repmat(obj.ph(:,1),1,obj.Nf);
                cosEl = sqrt(1 - sin(TH).^2.*sin(PH).^2);
                obj.E1 = (1./cosEl).*(cos(PH).*Eth - cos(TH).*sin(PH).*Eph);
                obj.E2 = (1./cosEl).*(cos(TH).*sin(PH).*Eth + cos(PH).*Eph);
                obj.E3 = zeros(size(obj.E1));
                obj.polBase = 'Ludwig2AE';
                obj = setEnames(obj);
                % Sort out singularities poles
                [~,iPole] = ismember([deg2rad(90),deg2rad(90);deg2rad(90),deg2rad(270)],[obj.th,obj.ph],'rows');
                iPole = iPole(iPole>0);
                [obj.E1(iPole,:),obj.E2(iPole,:),obj.E3(iPole,:)] = deal(0);
            end
        end
        
        function obj = pol2Ludwig2EA(obj)
            if ~strcmp(obj.polType,'Ludwig2EA')
                [Eth, Eph] = getEthEph(obj);
                TH = repmat(obj.th(:,1),1,obj.Nf);
                PH = repmat(obj.ph(:,1),1,obj.Nf);
                cosAl = sqrt(1 - sin(TH).^2.*cos(PH).^2);
                obj.E1 = (1./cosAl).*(cos(TH).*cos(PH).*Eth - sin(PH).*Eph);
                obj.E2 = (1./cosAl).*(sin(PH).*Eth + cos(TH).*cos(PH).*Eph);
                obj.E3 = zeros(size(obj.E1));
                obj.polBase = 'Ludwig2EA';
                obj = setEnames(obj);
                % Sort out singularities poles
%                 keyboard;
                [~,iPole] = ismember([deg2rad(90),deg2rad(0);deg2rad(90),deg2rad(180);deg2rad(90),deg2rad(360)],[obj.th,obj.ph],'rows');
                iPole = iPole(iPole>0);
                [obj.E1(iPole,:),obj.E2(iPole,:),obj.E3(iPole,:)] = deal(0);
            end
        end
        
        function obj = pol2Ludwig3(obj)
            if ~strcmp(obj.polType,'Ludwig3')
                [Eth, Eph] = getEthEph(obj);
                PH = repmat(obj.ph(:,1),1,obj.Nf);
                obj.E1 = cos(PH).*Eth - sin(PH).*Eph;
                obj.E2 = sin(PH).*Eth + cos(PH).*Eph;
                obj.E3 = zeros(size(obj.E1));
                obj.polBase = 'Ludwig3';
                obj = setEnames(obj);
            end
        end
        
        % Now change of polType - same idea as above with polBase
        % First get the linear pol and then go to what's required
        function [E1lin, E2lin] = getElin(obj)
            switch obj.polType
                case 'linear'
                    E1lin = obj.E1;
                    E2lin = obj.E2;
                case 'circular'
                    Del = 2*1i;
                    E1lin = sqrt(2)./Del.*(1i.*obj.E1 + 1i.*obj.E2);
                    E2lin = sqrt(2)./Del.*(-obj.E1 + obj.E2);
                case 'slant'
                    PSI = ones(size(obj.E1)).*obj.slant;
                    Del = 1;
                    E1lin = 1./Del.*(cos(PSI).*obj.E1 + sin(PSI).*obj.E2);
                    E2lin = 1./Del.*(-sin(PSI).*obj.E1 + cos(PSI).*obj.E2);
                otherwise
                    error(['Unknown polType property: ', obj.polType])
            end
        end
        
        function obj = pol2linear(obj)
            if ~strcmp(obj.polType,'linear')
                [E1lin, E2lin] = getElin(obj);
                obj.E1 = E1lin;
                obj.E2 = E2lin;
                obj.E3 = zeros(size(obj.E1));
                obj.polType = 'linear';
                obj = setEnames(obj);
            end
        end
        
        function obj = pol2circular(obj)
            if ~strcmp(obj.polType,'circular')
                [E1lin, E2lin] = getElin(obj);
                obj.E1 = 1/sqrt(2).*(E1lin - 1i.*E2lin);
                obj.E2 = 1/sqrt(2).*(E1lin + 1i.*E2lin);
                obj.polType = 'circular';
                obj = setEnames(obj);
            end
        end
        
        function obj = pol2slant(obj)
            if ~strcmp(obj.polType,'slant')
                [E1lin, E2lin] = getElin(obj);
                PSI = ones(size(obj.E1)).*obj.slant;
                obj.E1 = cos(PSI).*E1lin - sin(PSI).*E2lin;
                obj.E2 = sin(PSI).*E1lin + cos(PSI).*E2lin;
                obj.polType = 'slant';
                obj = setEnames(obj);
            end
        end
        
        %% Plotting methods
        
        plot(FF,varargin)
        
    end
    
    % Internal helper functions
    methods (Access = private)
        
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
                    switch obj.polBase
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
                        otherwise
                            error(['Unknown polBase property: ', obj.polBase]);
                    end
                    
                otherwise
                    error(['Unknown polType property: ', obj.polType]);
            end
        end
        
    end
    
end