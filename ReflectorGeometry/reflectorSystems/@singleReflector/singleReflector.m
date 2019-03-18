classdef singleReflector
    % Exactly the same as the GRASP single reflector wizard
    properties
        D(1,1) double {mustBeReal, mustBeFinite} = 10 % Reflector diameter in (m)
        FoD(1,1) double {mustBeReal, mustBeFinite} = 0.5 % Reflector F/D
        hoD(1,1) double {mustBeReal, mustBeFinite} = 0 % Reflector Offset/D
    end
    
    properties (SetAccess = private)
        F             % Reflector focal length
        h             % Reflector offset
        Dp            % Clearance between inner edge and the z-axis
        th_0          % Cone tilt angle from negative z-axis
        th_e          % Cone half angle (th* in GRASP technical description)
        th_f          % Feed tilt angle - points at middle of the reflector
        th_c          % Reflector rim angle with z-axis
        PR_chordX     % chord length along x-axis (a in GRASP technical description)
        PR_chordY     % chord length along y-axis
        apArea        % Projected aperture area
        P2            % Point at the edge of the PR on the x-axis
        P1            % Point at the edge of the PR on the x-axis
        P0            % Point at the centre of the PR (feed points here)
        PR            % Primary reflector object
        feedCoor      % Feed coordinate system
        apCoor        % Aperture coordinate system - located at centre of PR
    end
    
    methods
        function obj = singleReflector(D,FoD,hoD)
            if nargin > 0, obj.D = D; end
            if nargin > 1, obj.FoD = FoD; end
            if nargin > 2, obj.hoD = hoD; end
            
            % Basic geometry
            obj.F = obj.FoD*obj.D;
            obj.h = obj.hoD.*obj.D;
            obj.Dp = obj.h - obj.D/2;
            obj.th_0 = atan((2.*obj.F.*(obj.D + 2.*obj.Dp))./(4.*obj.F.^2 - obj.Dp.*(obj.D + obj.Dp)));
            obj.th_e = atan(2.*obj.F.*obj.D./(4.*obj.F.^2 + obj.Dp.*(obj.D + obj.Dp)));
            obj.th_f = 2.*atan((obj.Dp + obj.D/2)./(2.*obj.F));
            obj.th_c = atan(2.*obj.F./(obj.Dp + obj.D/2));
            obj.PR_chordX = obj.D./sin(obj.th_c);
            obj.PR_chordY = obj.D;
            obj.apArea = pi*(obj.D/2)^2;
            
            % Make the reflector
            surface = paraboloid(pnt3D(0,0,0),obj.F);
            rim = ellipticalRim([obj.h;0],[obj.D/2;obj.D/2]);
            PRcoor = coordinateSystem(pnt3D(0,0,0));
            obj.PR = reflector(surface,rim,PRcoor);
            
            % Calculate the extreme points of the reflector
            xy_xMin = [obj.Dp;0];
            xy_xMax = [obj.Dp+obj.D;0];
            xy_0 = [obj.h;0];
            xy = [xy_xMin,xy_xMax,xy_0];
            zEx = obj.PR.surface.getZ(xy(1,:),xy(2,:));
            obj.P1 = pnt3D(xy_xMin(1),xy_xMin(2),zEx(1));
            obj.P2 = pnt3D(xy_xMax(1),xy_xMax(2),zEx(2));
            obj.P0 = pnt3D(xy_0(1),xy_0(1),zEx(3));
            
            % Define the feed and aperture coordinates
            obj.feedCoor = coordinateSystem(pnt3D(0,0,obj.F));
            obj.feedCoor = obj.feedCoor.rotGRASP([pi-obj.th_f,0,pi]);
            obj.apCoor = coordinateSystem(pnt3D(obj.h,0,obj.P2.z + obj.F./2));
        end
        
        function [rho,drho_dth] = getThRhoMapping(obj,th)
            % Returns the th->rho mapping and its derivative 
            rho = 2.*obj.F.*tan(th./2);
            drho_dth = obj.F.*sec(th./2).^2;
        end
        
        function [pAp,pRefl] = rayTrace(obj,ph_in,th_in)
            % Returns the aperture and reflector points for a general 
            % element ray trace 
            [pRefl,reflectDir] = obj.PR.reflectRays(obj.feedCoor,ph_in,th_in);
            % Get the intersection of the reflection and the aperture plane
            % Provide the feed position as point in aperture plane
            n = repmat([0;0;1],1,size(reflectDir,2));    % Normal vector on aperture plane
            % Distance from intersection points to plane
            d = dot((obj.apCoor.origin.pointMatrix - pRefl.pointMatrix),n)./(dot(reflectDir,n));
            reflectDirScale = bsxfun(@times,reflectDir,d);
            pAp = pRefl.addVect(reflectDirScale);
        end
        
        function pathLengthStruct = getPathLength(obj,ph_in,th_in)
            % Calculate path length from:
            % feed to PR: FP
            % from PR to aperture: PA
            % from feed to aperture: FA
            if nargin == 1
                th_in = linspace(-obj.th_e,obj.th_e,21);
                ph_in = zeros(size(th_in));
            end
            ph_in = ph_in(:).';
            th_in = th_in(:).';
            [pAp,pRefl] = rayTrace(obj,ph_in,th_in);
            delFR = pRefl - obj.PR.coor.origin;
            delPA = pAp - pRefl;
            pathLengthStruct.ph = ph_in;
            pathLengthStruct.th = th_in;
            pathLengthStruct.FP = delFR.r;
            pathLengthStruct.PA = delPA.r;
            pathLengthStruct.FA = pathLengthStruct.FP + pathLengthStruct.PA;
        end
        
        function [FFM_F,MaskPointing,M] = getMask(obj,A)
            % Returns the reflector mask, from feed to PR, as a FarField
            % object. Also returns a matrix of pointing directions, as 
            % coordinate system objects, of the final rays after reflection 
            % throught the whole system. These can be used to assign
            % background temperature whan calculating antenna temperature.
            % Finally, the logical Mask is returned - true if being masked.
            % A can be a matrix of [ph,th] pairs, or a FarField object.  If
            % it is a FarField object it will be converted to a PhTh grid,
            % and those angles will be used
            if nargin == 1
                A = FarField;
            end
            if isa(A,'FarField')
                freq = A.freq;
            else
                freq = 1;
            end
            
            [M,ph_in,th_in] = obj.PR.getMask(obj.feedCoor,A);
            P = repmat(double(M(:)),1,numel(freq));
            FFM_F = FarField.farFieldFromPowerPattern(ph_in(:),th_in(:),P,freq);
            % Build the pointing matrix
            MaskPointing(size(ph_in)) = coordinateSystem;
            % First those outside mask - centered at feed currently
            [MaskPointing(~M).base] = deal(obj.feedCoor);
            % Those in mask - origin at global base and pointing up so do
            % nothing...
            % Fix the non-masked position pointing angles - in the global
            % coordinate system
            maskI = find(~M);
            for mm = maskI
                MaskPointing(mm) = MaskPointing(mm).rotGRASP([th_in(mm),ph_in(mm),0]);
                MaskPointing(mm) = MaskPointing(mm).getInGlobal;    % Rotate to global Coor
                MaskPointing(mm).origin = pnt3D;    % Force to centre of global Coor
            end
        end
        
        %% Plotting
        function plot(obj,N)
            if nargin == 1
                N = 101;
            end
            lineWidthRefl = 2;
            lineWidthRays = 0.5;
            % Plot the reflector surface in the symmetry plane
            [surfPoints] = obj.PR.getPointCloud(N,'x0');
            x = surfPoints.x;
            z = surfPoints.z;
            plot(x,z,'k','linewidth',lineWidthRefl), hold on, grid on
            % Plot the feed point
            plot(obj.feedCoor.origin.x,obj.feedCoor.origin.z,'k.','markersize',1.5*obj.PR_chordX)
            % Plot the edge rays
            xR = [obj.P1.x,obj.P1.x,obj.feedCoor.origin.x,obj.P2.x,obj.P2.x];
            zR = [obj.apCoor.origin.z,obj.P1.z,obj.feedCoor.origin.z,obj.P2.z,obj.apCoor.origin.z];
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
            obj.PR.plot(N)
            if coorFlag(1), obj.feedCoor.plot(obj.D/10); end
            if coorFlag(2), obj.apCoor.plot(obj.D/10); end
            if coorFlag(3), obj.PR.coor.plot(obj.D/10); end
            xlabel('x-axis (m)')
            ylabel('y-axis (m)')
            zlabel('z-axis (m)')
        end
        
        function plotRayTrace(obj,Nray,Nrefl)
            if nargin == 1
                Nray = 21;
                Nrefl = 101;
            elseif nargin == 2
                Nrefl = 101;
            end
            th_in = linspace(-obj.th_e,obj.th_e,Nray);
            ph_in = zeros(size(th_in));
            obj.plot3D(Nrefl)
            hold on
            [pAp,pRefl] = rayTrace(obj,ph_in,th_in);
            rayColor = ones(1,3).*0;
            rayWidth = 0.5;
            plotLines(obj.feedCoor.origin,pRefl,'lineColor',rayColor,'lineWidth',rayWidth)
            plotLines(pRefl,pAp,'lineColor',rayColor,'lineWidth',rayWidth)
            view([0,0])
        end
        
    end
    
end