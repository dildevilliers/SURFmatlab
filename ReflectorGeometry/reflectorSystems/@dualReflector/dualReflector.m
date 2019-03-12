classdef dualReflector
    % Exactly the same as the GRASP dual reflector wizard
    % Use the variable names of Granet 2002 paper
    properties
        Dm(1,1) double {mustBeReal, mustBeFinite} = 10 
        Lm(1,1) double {mustBeReal, mustBeFinite} = 1  
        Ls(1,1) double {mustBeReal, mustBeFinite} = 2.5  
        th_e(1,1) double {mustBeReal, mustBeFinite} = deg2rad(15) 
        th_0(1,1) double {mustBeReal, mustBeFinite} = 0
        beta(1,1) double {mustBeReal, mustBeFinite} = 0 
        sigma(1,1) double {mustBeReal, mustBeFinite} = 1 
        th_ext(1,1) double {mustBeReal, mustBeFinite} = 0 % Extension angle in (rad)
        symFact_ext(1,1) double {mustBeReal, mustBeFinite} = 0 % Symmetry factor of the SR extension. 0 is symmetric, 1 is bottom, and -1 is top extension
    end
    
    properties (SetAccess = private)
        % Granet 2002 paper for description
        % PR is used instead of MR
        type % String describing the system
        F
        a
        f
        e
        alpha
        h
        th_U
        th_L
        Dsx
        Dsy
        dSRPR
        dFPR
        Ht
        Lt
        C_SR
        Feq     % Equivalent focal length
        % Extreme points on the reflectors and aperture
        P0
        P1
        P2
        Q0
        Q1
        Q2
        R0
        R1
        R2
        % Physical information
        PR_chordX
        PR_chordY
        SR_chordX
        SR_chordY
        apArea        % Projected aperture area
        % Actual objects
        PR  % Primary reflector
        SR  % Secondary reflector
        feedCoor
        apCoor
    end
    
    methods
        % Constructor is over specified - so problems can occur when called
        % directly.  This is to maintain flexibility with defining
        % symmetric and offset systems in one class.  Use one of the
        % meta-constructors below, which are specified according to the
        % Granet paper options and calls this constructor correctly
        function obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext)
            if nargin >= 3 % Set 6, symmetric, min block
                obj.Dm = Dm;
                obj.Lm = Lm;
                obj.th_e = th_e;
            end
            if nargin >= 4 % Set 1, symmetric, no block condition
                obj.Ls = Ls;
            end
            if nargin >= 5 % Set 8 offset
                obj.th_0 = th_0;
                obj.beta = beta;
            end
            if nargin >= 7
                obj.sigma = sigma;
            end
            if nargin >= 8, obj.th_ext = th_ext; end
            if nargin == 9, obj.symFact_ext = symFact_ext; end
            
            if obj.sigma == 1
                obj.type = ['Gregorian'];
                SRhandle = str2func('ellipsoid');
            elseif obj.sigma == -1
                obj.type = ['Cassegrain'];
                SRhandle = str2func('hyperboloid');
            end

            % Build the geometry
            if obj.th_0 == 0 && obj.beta == 0
                % Symmetrical case
                obj.type = ['Symmetric ',obj.type];
                obj.f = obj.Ls*(obj.sigma.*obj.Dm - 4.*obj.Lm.*tan(obj.th_e/2))./(2.*obj.sigma.*obj.Dm + 8.*obj.Ls.*tan(obj.th_e/2));   % (4)
                obj.F = obj.Lm + 2.*obj.f;  % (5)
                obj.a = obj.Ls - obj.f; % (6)
                obj.Dsx = 4.*(obj.Ls - obj.f)./(1./sin(obj.th_e) + obj.sigma.*(16.*obj.F.^2 + obj.Dm.^2)./(8.*obj.F.*obj.Dm));   % (7)
                obj.e = obj.f./obj.a;
                obj.alpha = 0;
                obj.th_U = -2.*atan(obj.Dm/(4.*obj.F));
            else % Offset Case
                obj.type = ['Offset ',obj.type];
                obj.e = (1 - obj.sigma*sqrt(tan(obj.beta/2)/tan((obj.beta - obj.th_0)/2)))./(1 + obj.sigma*sqrt(tan(obj.beta/2)/tan((obj.beta - obj.th_0)/2)));  % (5)
                obj.alpha = 2*atan((obj.e+1)/(obj.e-1)*tan(obj.beta/2));   % (6)
                obj.a = obj.Ls/(2 + (obj.e^2-1)/(obj.e*cos(obj.beta-obj.th_0) + 1));  % (22)
                obj.f = obj.a*obj.e;    % (15)
                obj.th_U = 2*atan((1+obj.e)/(1-obj.e)*tan((obj.alpha-obj.sigma*obj.th_e)/2)) + obj.beta;   % (3)
                obj.F = obj.Dm/(4*(tan(-obj.th_U/2) - tan(-obj.th_0/2)));   % (26)
            end
            obj.h = 2*obj.F*tan(-obj.th_0/2);   % (23)
            obj.th_L = -2*atan((2*obj.h-obj.Dm)/(4*obj.F)); % (4)
            obj.Dsx = -obj.sigma*obj.a*(((obj.e^2-1)*sin(obj.beta-obj.th_U))/(obj.e*cos(obj.beta-obj.th_U)+1) - ((obj.e^2-1)*sin(obj.beta-obj.th_L))/(obj.e*cos(obj.beta-obj.th_L)+1));    % (27)
            obj.dSRPR = obj.h - obj.Dm/2 - obj.a*(obj.sigma-1)/2*(obj.e^2-1)*sin(obj.th_U)/(obj.e*cos(obj.beta-obj.th_U)+1) + obj.a*(obj.sigma+1)/2*(obj.e^2-1)*sin(obj.th_L)/(obj.e*cos(obj.beta-obj.th_L)+1); % (11)
            obj.dFPR = obj.h - obj.Dm/2 + 2*obj.f*sin(obj.beta);   % (10)
            obj.Lt = -obj.a*((obj.sigma + 1)/2)*(((obj.e^2-1)*cos(obj.th_L))/(obj.e*cos(obj.beta-obj.th_L) + 1)) + obj.a*((obj.sigma - 1)/2)*(((obj.e^2-1)*cos(obj.th_U))/(obj.e*cos(obj.beta-obj.th_U) + 1)) - (2*obj.h - obj.Dm)^2/(16*obj.F) + obj.F;  % (12)
            obj.Ht = obj.h + obj.Dm/2 - obj.a*((obj.sigma - 1)/2)*(((obj.e^2-1)*sin(obj.th_L))/(obj.e*cos(obj.beta-obj.th_L) + 1)) + obj.a*((obj.sigma + 1)/2)*(((obj.e^2-1)*sin(obj.th_U))/(obj.e*cos(obj.beta-obj.th_U) + 1));  % (13)
            if obj.th_0 == 0 && obj.beta == 0
                obj.Ht = obj.Dm;    % Special case for symmetric - Ht gets new meaning
            end
            FoD = cot(obj.th_e/2)/4;
            obj.Feq = FoD.*obj.Dm;
            
            % Feed coordinate system needed later and in global coordinates
            feedCoor = coordinateSystem();
            feedCoor.origin = pnt3D(-2*obj.f*sin(obj.beta),0,-2*obj.f*cos(obj.beta)); 
            obj.feedCoor = feedCoor.rotGRASP([obj.beta+obj.alpha,0,0]);

            % The important points
            obj.Q0 = pnt3D(obj.h,0,obj.h^2/(4*obj.F) - obj.F);
            obj.Q1 = pnt3D(obj.h-obj.Dm/2,0,(2*obj.h - obj.Dm)^2/(16*obj.F) - obj.F);
            obj.Q2 = pnt3D(obj.h+obj.Dm/2,0,(2*obj.h + obj.Dm)^2/(16*obj.F) - obj.F);
            
            Rz = obj.Dm/2;
            obj.R0 = obj.Q0;
            obj.R1 = obj.Q1;
            obj.R2 = obj.Q2;
            [obj.R0.z,obj.R1.z,obj.R2.z] = deal(Rz);
            
            OP0 = obj.sigma*(2*obj.a-obj.Ls);   % (32)
            OP1 = -obj.sigma*obj.a*(obj.e^2-1)/(obj.e*cos(obj.th_L - obj.beta) + 1); % (33)
            OP2 = -obj.sigma*obj.a*(obj.e^2-1)/(obj.e*cos(obj.th_U - obj.beta) + 1); % (34)
            
            obj.P0 = pnt3D(obj.sigma*OP0*sin(obj.th_0),0,obj.sigma*OP0*cos(obj.th_0));
            obj.P1 = pnt3D(obj.sigma*OP1*sin(obj.th_L),0,obj.sigma*OP1*cos(obj.th_L));
            obj.P2 = pnt3D(obj.sigma*OP2*sin(obj.th_U),0,obj.sigma*OP2*cos(obj.th_U));
            
            Cx = (distanceCart(obj.feedCoor.origin,obj.P1)*sin(obj.alpha + obj.sigma.*obj.th_e) + distanceCart(obj.feedCoor.origin,obj.P2).*sin(obj.alpha - obj.sigma.*obj.th_e))/2; % (38)
            Cy = 0; % (38)
            Cz = obj.a*sqrt(1 + Cx^2/(obj.f^2 - obj.a^2)) - obj.f; % (38)
            obj.C_SR = pnt3D(Cx,Cy,Cz);
            
            ph = linspace(0,2*pi,1001);
            obj.Dsy = max((2*obj.a*(obj.e^2-1)*sin(obj.th_e).*sin(ph))./(obj.e.*(-sin(obj.alpha).*sin(obj.th_e).*cos(ph) + cos(obj.alpha).*cos(obj.th_e)) - 1));   % (39)

            % Derived geometry
            obj.PR_chordX = distanceCart(obj.Q2,obj.Q1);
            obj.PR_chordY = obj.Dm;
            obj.SR_chordX = distanceCart(obj.P2,obj.P1);
            obj.SR_chordY = obj.Dsy;
            obj.apArea = pi*(obj.Dm/2)^2;

            % Build the reflectors and coordinate systems
            obj.SR = reflector;
            obj.SR.surface = SRhandle(2*obj.a,2*obj.f);
            obj.SR.rim = ellipticalRim([obj.C_SR.x,obj.C_SR.y],[obj.Dsx,obj.Dsy]./2);
            SRcoor = coordinateSystem;
            obj.SR.coor = SRcoor.rotGRASP([obj.beta,0,0]);
            
            obj.PR = reflector;
            obj.PR.surface = paraboloid(pnt3D(0,0,-obj.F),obj.F);
            obj.PR.rim = ellipticalRim([obj.Q0.x,obj.Q0.y],[obj.Dm,obj.Dm]./2);
            obj.PR.coor = coordinateSystem();
            
            obj.apCoor = coordinateSystem(obj.R0);
        end
        
        function [rho,drho_dth] = getThRhoMapping(obj,th)
            % Returns the th->rho mapping and its derivative 
            rho = 2.*obj.Feq.*tan(th./2);
            drho_dth = obj.Feq.*sec(th./2).^2;
        end
        
        function [pAp,pRefl] = rayTrace(obj,ph_in,th_in)
            % Returns the aperture and reflector points for an Nray 
            % element ray trace in the y=0 plane 
%             [pRefl,reflectDir] = obj.PR.reflectRays(obj.feedCoor,ph_in,th_in);
%             % Get the intersection of the reflection and the aperture plane
%             % Provide the feed position as point in aperture plane
%             n = repmat([0;0;1],1,size(reflectDir,2));    % Normal vector on aperture plane
%             % Distance from intersection points to plane
%             d = dot((obj.apCoor.origin.pointMatrix - pRefl.pointMatrix),n)./(dot(reflectDir,n));
%             reflectDirScale = bsxfun(@times,reflectDir,d);
%             pAp = pRefl.addVect(reflectDirScale);
        end
        
        function pathLengthStruct = getPathLength(obj,ph_in,th_in)
            % Calculate path length from:
            % feed to PR: FP
            % from PR to aperture: PA
            % from feed to aperture: FA
%             if nargin == 1
%                 th_in = linspace(-obj.th_e,obj.th_e,21);
%                 ph_in = zeros(size(th_in));
%             end
%             ph_in = ph_in(:).';
%             th_in = th_in(:).';
%             [pAp,pRefl] = rayTrace(obj,ph_in,th_in);
%             delFR = pRefl - obj.PR.coor.origin;
%             delPA = pAp - pRefl;
%             pathLengthStruct.ph = ph_in;
%             pathLengthStruct.th = th_in;
%             pathLengthStruct.FP = delFR.r;
%             pathLengthStruct.PA = delPA.r;
%             pathLengthStruct.FA = pathLengthStruct.FP + pathLengthStruct.PA;
        end
        
        function FFM_F = getMask(obj,A)
            % Returns the reflector mask, from feed to PR, as a FarField
            % object.  
            % A can be a matrix of [ph,th] pairs, or a FarField object.  If
            % it is a FarField object it will be converted to a PhTh grid,
            % and those angles will be used
%             if nargin == 1
%                 A = FarField;
%             end
%             if isa(A,'FarField')
%                 freq = A.freq;
%             else
%                 freq = 1;
%             end
%             
%             [M,ph_in,th_in] = getMask(obj.PR,obj.feedCoor,A);
%             P = repmat(double(M(:)),1,numel(freq));
%             FFM_F = FarField.farFieldFromPowerPattern(ph_in(:),th_in(:),P,freq);
%             FFM_F = FFM_F.setXrange('pos');
%             FFM_F = FFM_F.currentForm2Base;
        end
        
        %% Plotting
        function plot(obj,N)
            if nargin == 1
                N = 101;
            end
            lineWidthRefl = 2;
            lineWidthRays = 0.5;
            % Plot the reflector surface in the symmetry plane
            [surfPointsPR] = obj.PR.getPointCloud(N,'x0');
            [surfPointsSR] = obj.SR.getPointCloud(N,'x0');
            plot(surfPointsPR.x,surfPointsPR.z,'k','linewidth',lineWidthRefl), hold on, grid on
            plot(surfPointsSR.x,surfPointsSR.z,'k','linewidth',lineWidthRefl)
            % Plot the feed point
            plot(obj.feedCoor.origin.x,obj.feedCoor.origin.z,'k.','markersize',1.0*obj.SR_chordX)
            % Plot the edge rays
            xR = [obj.R2.x,obj.Q2.x,obj.P2.x,obj.feedCoor.origin.x,obj.P1.x,obj.Q1.x,obj.R1.x];
            zR = [obj.R2.z,obj.Q2.z,obj.P2.z,obj.feedCoor.origin.z,obj.P1.z,obj.Q1.z,obj.R1.z];
            plot(xR,zR,'k','linewidth',lineWidthRays)
            xlabel('x (m)')
            ylabel('z (m)')
            axis equal
        end
        
        function plot3D(obj,N,coorFlag)
            % N is the approximated root of the number of points
            % coorFlag is a vector of length 3 indicating if the coordinate
            % systems should be plotted: [feedCoor,apertureCoor,reflectorCoor]
            % Default is [1,0,0]
            if nargin == 1
                N = 10000;
                coorFlag = [1,0,0];
            elseif nargin  == 2
                coorFlag = [1,0,0];
            end
            obj.SR.plot(N), hold on
            obj.PR.plot(N)
            if coorFlag(1), obj.feedCoor.plot(obj.Dm/10); end
            if coorFlag(2), obj.apCoor.plot(obj.Dm/10); end
            if coorFlag(3)
                obj.PR.coor.plot(obj.Dm/10); 
                obj.SR.coor.plot(obj.Dm/10);
            end
            xlabel('x-axis (m)')
            ylabel('y-axis (m)')
            zlabel('z-axis (m)')
        end
        
        function plotRayTrace(obj,Nray,Nrefl)
%             if nargin == 1
%                 Nray = 21;
%                 Nrefl = 101;
%             elseif nargin == 2
%                 Nrefl = 101;
%             end
%             th_in = linspace(-obj.th_e,obj.th_e,Nray);
%             ph_in = zeros(size(th_in));
%             obj.plot(Nrefl)
%             hold on
%             [pAp,pRefl] = rayTrace(obj,ph_in,th_in);
%             rayColor = ones(1,3).*0;
%             rayWidth = 0.5;
%             plotLines(obj.feedCoor.origin,pRefl,'lineColor',rayColor,'lineWidth',rayWidth)
%             plotLines(pRefl,pAp,'lineColor',rayColor,'lineWidth',rayWidth)
%             view([0,0])
        end
        
    end
    
end