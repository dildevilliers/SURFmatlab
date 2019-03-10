classdef symmetricParaboloid
    properties
        D(1,1) double {mustBeReal, mustBeFinite} = 10 % in (m)
        FoD(1,1) double {mustBeReal, mustBeFinite} = 0.5 % F/D
        feedCoor(1,1) coordinateSystem = coordinateSystem(pnt3D(0,0,0),[1;0;0],[0;-1;0]);
        apertureCoor(1,1) coordinateSystem = coordinateSystem();
    end
    
    properties (SetAccess = private)
        F
        th_e
        P_PR_xMax     % Point at the edge of the PR on the x-axis
        P_PR_xMin     % Point at the edge of the PR on the x-axis
        P_PR_yMax     % Point at the edge of the PR on the y-axis
        P_PR_yMin     % Point at the edge of the PR on the y-axis
        P_PR0         % Point at the vertex of the PR (feed points here)
        chordX        % chord length along x-axis
        chordY        % chord length along y-axis
        PR            % Primary reflector object
        apArea        % Projected aperture area
    end
    
    methods
        function obj = symmetricParaboloid(D,FoD,feedCoor,apertureCoor)
            if nargin > 0, obj.D = D; end
            if nargin > 1, obj.FoD = FoD; end
            if nargin > 2, obj.feedCoor = feedCoor; end
            if nargin > 3, obj.apertureCoor = apertureCoor; end
            
            obj.F = obj.FoD*obj.D;
            surface = paraboloid(pnt3D(0,0,-obj.F),obj.F);
            rim = ellipticalRim([0;0],[obj.D/2;obj.D/2]);
            PRcoor = coordinateSystem(pnt3D(0,0,0));
            obj.PR = reflector(surface,rim,PRcoor);
            obj.th_e = fd2th0(obj.FoD);
            % Calculate the extreme points of the reflector
            xy_xMin = [-obj.D/2;0];
            xy_xMax = [obj.D/2;0];
            xy_yMin = [0;-obj.D/2];
            xy_yMax = [0;obj.D/2];
            xy_0 = [0;0];
            xy = [xy_xMin,xy_xMax,xy_yMin,xy_yMax,xy_0];
            zEx = obj.PR.surface.getZ(xy(1,:),xy(2,:));
            obj.P_PR_xMin = pnt3D(xy_xMin(1),xy_xMin(2),zEx(1));
            obj.P_PR_xMax = pnt3D(xy_xMax(1),xy_xMax(2),zEx(2));
            obj.P_PR_yMin = pnt3D(xy_yMin(1),xy_yMin(2),zEx(3));
            obj.P_PR_yMax = pnt3D(xy_yMax(1),xy_yMax(2),zEx(4));
            obj.P_PR0 = pnt3D(xy_0(1),xy_0(1),zEx(5));
            % Derived physical quantities
            obj.chordX = obj.D;
            obj.chordY = obj.D;
            obj.apArea = pi*(obj.D/2)^2;
        end
        
        function [rho,drho_dth] = getThRhoMapping(obj,th)
            % Returns the th->rho mapping and its derivative 
            rho = 2.*obj.F.*tan(th./2);
            drho_dth = obj.F.*sec(th./2).^2;
        end
        
        function [pAp,pRefl] = rayTrace(obj,ph_in,th_in)
            % Returns the aperture and reflector points for an Nray 
            % element ray trace in the y=0 plane 
            [pRefl,reflectDir] = obj.PR.reflectRays(obj.feedCoor,ph_in,th_in);
            % Get the intersection of the reflection and the aperture plane
            % Provide the feed position as point in aperture plane
            n = repmat([0;0;1],1,size(reflectDir,2));    % Normal vector on aperture plane
            % Distance from intersection points to plane
            d = dot((obj.feedCoor.origin.pointMatrix - pRefl.pointMatrix),n)./(dot(reflectDir,n));
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
        
        function FFM_F = getMask(obj,A)
            % Returns the reflector mask, from feed to PR, as a FarField
            % object.  
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
            
            [M,ph_in,th_in] = getMask(obj.PR,obj.feedCoor,A);
            P = repmat(double(M(:)),1,numel(freq));
            FFM_F = FarField.farFieldFromPowerPattern(ph_in(:),th_in(:),P,freq);
            FFM_F = FFM_F.setXrange('pos');
            FFM_F = FFM_F.currentForm2Base;
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
            y = surfPoints.y;
            z = surfPoints.z;
            plot3(x,y,z,'k','linewidth',lineWidthRefl), hold on, grid on
            % Plot the feed point
            plot3(obj.feedCoor.origin.x,obj.feedCoor.origin.y,obj.feedCoor.origin.z,'k.','markersize',1.5*obj.chordX)
            % Plot the edge rays
            xR = [obj.P_PR_xMin.x,obj.P_PR_xMin.x,obj.feedCoor.origin.x,obj.feedCoor.origin.x,obj.feedCoor.origin.x,obj.P_PR_xMax.x,obj.P_PR_xMax.x];
            yR = zeros(size(xR));
            zR = [obj.apertureCoor.origin.z,obj.P_PR_xMin.z,obj.feedCoor.origin.z,obj.P_PR0.z,obj.feedCoor.origin.z,obj.P_PR_xMax.z,obj.apertureCoor.origin.z];
            plot3(xR,yR,zR,'k','linewidth',lineWidthRays)
            xlabel('x (m)')
            zlabel('z (m)')
            axis equal
            view([0,0])
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
            if coorFlag(2), obj.apertureCoor.plot(obj.D/10); end
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
            obj.plot(Nrefl)
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