classdef symmetricParaboloid
    properties
        D(1,1) double {mustBeReal, mustBeFinite} = 10 % in (m)
        FoD(1,1) double {mustBeReal, mustBeFinite} = 0.5 % F/D
        feedCoor(1,1) coordinateSystem = coordinateSystem([0;0;0],[1;0;0],[0;-1;0]);
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
            surface = paraboloid([0;0;-obj.F],obj.F);
            rim = ellipticalRim([0;0],[obj.D/2;obj.D/2]);
            PRcoor = coordinateSystem([0,0,0]);
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
            obj.P_PR_xMin = [xy_xMin;zEx(1)];
            obj.P_PR_xMax = [xy_xMax;zEx(2)];
            obj.P_PR_yMin = [xy_yMin;zEx(3)];
            obj.P_PR_yMax = [xy_yMax;zEx(4)];
            obj.P_PR0 = [xy_0;zEx(5)];
            % Derived physical quantities
            obj.chordX = norm(obj.P_PR_xMax - obj.P_PR_xMin);
            obj.chordY = norm(obj.P_PR_yMax - obj.P_PR_yMin);
            obj.chordX = obj.D;
            obj.chordY = obj.D;
            obj.apArea = pi*(obj.D/2)^2;
        end
        
        function [rho,drho_dth] = getThRhoMapping(obj,th)
            % Returns the th->rho mapping and its derivative 
            rho = 2.*obj.F.*tan(th./2);
            drho_dth = obj.F.*sec(th./2).^2;
        end
        
        function pathLengthStruct = getPathLength(obj,th,ph)
            % Calculate path length from:
            % feed to PR: FP
            % from PR to aperture: PA
            % from feed to aperture: FA
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
            x = surfPoints(1,:);
            z = surfPoints(3,:);
            plot(x,z,'k','linewidth',lineWidthRefl), hold on, grid on
            % Plot the feed point
            plot(obj.feedCoor.origin(1),obj.feedCoor.origin(3),'k.','markersize',1.5*obj.chordX)
            % Plot the edge rays
            xR = [obj.P_PR_xMin(1),obj.P_PR_xMin(1),obj.feedCoor.origin(1),obj.feedCoor.origin(1),obj.feedCoor.origin(1),obj.P_PR_xMax(1),obj.P_PR_xMax(1)];
            zR = [obj.apertureCoor.origin(3),obj.P_PR_xMin(3),obj.feedCoor.origin(3),obj.P_PR0(3),obj.feedCoor.origin(3),obj.P_PR_xMax(3),obj.apertureCoor.origin(3)];
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
                N = 101;
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
        
    end
    
end