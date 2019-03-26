classdef offsetParaboloid
    properties
        D(1,1) double {mustBeReal, mustBeFinite} = 10 % in (m)
        FoD(1,1) double {mustBeReal, mustBeFinite} = 0.5 % F/D
        feedCoor(1,1) coordinateSystem = coordinateSystem([0;0;0],[1;0;0],[0;-1;0]);
        apertureCoor(1,1) coordinateSystem = coordinateSystem();
        
        % Offset paraboloid object start list.
        MR;
        coor_Feed;
        % Offset paraboloid object end list.
    end
    
    properties (SetAccess = private)
        %F
        th_e
        P_PR_xMax     % Point at the edge of the PR on the x-axis
        P_PR_xMin     % Point at the edge of the PR on the x-axis
        P_PR_yMax     % Point at the edge of the PR on the y-axis
        P_PR_yMin     % Point at the edge of the PR on the y-axis
        P_PR0         % Point at the vertex of the PR (feed points here)
        chordX        % chord length along x-axis
        chordY        % chord length along y-axis
        PR            % Primary reflector object
        
        % Offset paraboloid variables start list. 
        F;
        alpha;
        beta;
        
        ReflectorDiameter;
        ReflectorFoD;
        ReflectorOffset;
        % Offset paraboloid variables end list.
    end
    
    methods
        function obj = offsetParaboloid(argFoD,argReflectorDiameter,argReflectorOffset) % SymParaArgs: D,FoD,feedCoor,apertureCoor
            if nargin > 0, obj.ReflectorFoD = argFoD; end
            if nargin > 1, obj.ReflectorDiameter = argReflectorDiameter; end
            if nargin > 2, obj.ReflectorOffset = argReflectorOffset; end
            %if nargin > 3, obj.apertureCoor = apertureCoor; end
            
            % Offset paraboloid alpha,beta and F assignment. The 3 unique
            % describing variables.
            
            % Offset paraboloid assignment ended.
            
            %obj.F = obj.FoD*obj.D;
            %surface = paraboloid([0;0;0],obj.F);
            %rim = ellipticalRim([0;0],[obj.D/2;obj.D/2]);
            %PRcoor = coordinateSystem([0,0,-obj.F]);
            %obj.PR = reflector(surface,rim,PRcoor);
            %obj.th_e = fd2th0(obj.FoD);
            %obj.P_PR_xMax = [obj.D/2;0;0];
            %obj.P_PR_xMin = [-obj.D/2;0;0];
            %obj.P_PR_yMax = [0;obj.D/2;0];
            %obj.P_PR_yMin = [0;-obj.D/2;0];
            %obj.chordX = norm(obj.P_PR_xMax - obj.P_PR_xMin);
            %obj.chordY = norm(obj.P_PR_yMax - obj.P_PR_yMin);
            %obj.P_PR0 = PRcoor.origin;
            
            % Offset Paraboloid PROPERTY ASSIGNMENT AREA START.
            obj.F = argReflectorDiameter*argFoD;
            rim_move_X_ = argReflectorDiameter*argReflectorOffset; %Alternate name: Xcentre.
            rim_radius_X_ = argReflectorDiameter/2; %Alternate name: XHalfAxis.
            rim_radius_Y_ = argReflectorDiameter/2; %Alternate name: YHalfAxis.
            % Offset Paraboloid PROPERTY ASSIGNMENT AREA END..
            
            MRParaboloid = paraboloid([0;0;0],obj.F);
            MRRim = ellipticalRim([rim_move_X_;0],[rim_radius_X_;rim_radius_Y_],0);
            MRCoor = coordinateSystem([0,0,0]);
            
            obj.MR = reflector(MRParaboloid,MRRim,MRCoor);
            
            %NOW CALCULATE ALPHA BETA ANGLES FROM GRASP INFORMATION BELOW:
            xDist_ = rim_move_X_;
            zCheck = MRParaboloid.getZ(rim_move_X_, 0);
            yDist_ = obj.F - zCheck;
            betaAngle = atan(xDist_/yDist_);
            obj.beta = betaAngle;
            
            Px = rim_move_X_ - rim_radius_X_;
            Pz = MRParaboloid.getZ(Px, 0);
            gammaAngle = atan(Px/(obj.F-Pz));
            alpha1 = betaAngle - gammaAngle;
            
            PxMAX = rim_move_X_ + rim_radius_X_;
            PzMAX = MRParaboloid.getZ(PxMAX, 0);
            omegaAngle = atan(PxMAX/(obj.F-PzMAX));
            alpha2 = omegaAngle - betaAngle;
            
            obj.alpha = (alpha1+alpha2)/2;
            %NOW CALCULATE ALPHA BETA ANGLES FROM GRASP INFORMATION ABOVE.
            
            % Offset Paraboloid OBJECT CONSTUCTION AREA START.
            
            obj.coor_Feed = coordinateSystem([0,0,obj.F],[1;0;0],[0;1;0],MRCoor); %
            obj.coor_Feed = obj.coor_Feed.rotGRASP([((180.0-obj.beta/pi*180.0)*pi/180),0,0]); %
            % Offset Paraboloid OBJECT CONSTUCTION AREA END.
        end
        
        function obj = calculateAlphaBetaAngles(obj)
            xDist_ = rim_move_X_;
            zCheck = obj.surface.getZ(rim_move_X_, 0);
            yDist_ = obj.F - zCheck;
            betaAngle = atan(yDist_/xDist_);
            obj.beta = betaAngle;
            
            Px = rim_move_X_ - rim_radius_X_;
            Pz = obj.surface.getZ(Px, 0);
            
            gammaAngle = atan(Px/(obj.F-Pz));
            obj.alpha = betaAngle - gammaAngle;
        end
        
        function obj = buildReflectorUsingAlphaBetaAngles(obj) %, argF, argAlpha, argBeta
            %if nargin > 0, obj.F = argF; end
            %if nargin > 1, obj.alpha = argAlpha; end
            %if nargin > 2, obj.beta = argBeta; end
            %if nargin > 3, obj.apertureCoor = apertureCoor; end
            
            % Offset paraboloid alpha,beta and F assignment. The 3 unique
            % describing variables.
            
            % Offset paraboloid assignment ended.
            
            %obj.F = obj.FoD*obj.D;
            %surface = paraboloid([0;0;0],obj.F);
            %rim = ellipticalRim([0;0],[obj.D/2;obj.D/2]);
            %PRcoor = coordinateSystem([0,0,-obj.F]);
            %obj.PR = reflector(surface,rim,PRcoor);
            %obj.th_e = fd2th0(obj.FoD);
            %obj.P_PR_xMax = [obj.D/2;0;0];
            %obj.P_PR_xMin = [-obj.D/2;0;0];
            %obj.P_PR_yMax = [0;obj.D/2;0];
            %obj.P_PR_yMin = [0;-obj.D/2;0];
            %obj.chordX = norm(obj.P_PR_xMax - obj.P_PR_xMin);
            %obj.chordY = norm(obj.P_PR_yMax - obj.P_PR_yMin);
            %obj.P_PR0 = PRcoor.origin;
            
            % Offset Paraboloid PROPERTY ASSIGNMENT AREA START.
            %pi = 3.141592654;
            Eccentricity = 1.001;
            
            a__ = ((cos(obj.beta)+Eccentricity*cos(obj.alpha))/sin(obj.beta))^2.0;
            a__SQRT_ = ((cos(obj.beta)+Eccentricity*cos(obj.alpha))/sin(obj.beta));
            asquared_ = (((1.0+Eccentricity)^2.0*(sin(obj.alpha))^2.0)/((1.0-Eccentricity^2.0)*(sin(obj.beta))^2.0+(cos(obj.beta)+Eccentricity*cos(obj.alpha))^2.0))*(obj.F^2.0);
            b__ = ((cos(obj.beta)+Eccentricity*cos(obj.alpha))/sin(obj.beta))*(((1.0+Eccentricity)*cos(obj.alpha))/sin(obj.beta));
            bsquared_ = (((1.0+Eccentricity)^2.0*(cos(obj.beta)+Eccentricity*cos(obj.alpha))^2.0*(sin(obj.alpha))^2.0)/(((1.0-Eccentricity^2.0)*(sin(obj.beta))^2.0+(cos(obj.beta)+Eccentricity*cos(obj.alpha))^2.0)^2.0))*(obj.F^2.0);
            c__ = (((1.0+Eccentricity)*cos(obj.alpha))/sin(obj.beta))^2.0;
            c__SQRT_ = (((1.0+Eccentricity)*cos(obj.alpha))/sin(obj.beta));

            L__ = (((1.0-Eccentricity^2.0)/((1.0-Eccentricity)^2.0))*obj.F^2.0);
            k__ = L__/((obj.F^2.0)/((1.0-Eccentricity)^2.0));
            M__ = (a__+k__);
            P__ = (obj.F^2.0*c__-L__+k__*(((Eccentricity/(1.0-Eccentricity))^2.0)*(obj.F^2.0)));
            R__ = (2.0*b__*obj.F-2.0*k__*(Eccentricity/(1.0-Eccentricity))*obj.F);
            
            ZQsolution1_ = ((-R__+sqrt(((R__)^2.0-4.0*M__*P__)))/(2.0*M__));
            ZQsolution2_ = ((-R__-sqrt(((R__)^2.0-4.0*M__*P__)))/(2.0*M__));
            YQsolution1_ = a__SQRT_*ZQsolution1_+c__SQRT_*obj.F;
            YQsolution2_ = a__SQRT_*ZQsolution2_+c__SQRT_*obj.F;
            y0_ = (((1.0+Eccentricity)*(Eccentricity*cos(obj.beta)+cos(obj.alpha))*sin(obj.beta))/((1.0-Eccentricity^2.0)*(sin(obj.beta))^2.0+(cos(obj.beta)+Eccentricity*cos(obj.alpha))^2.0))*obj.F;
            
            rim_move_X_ = (YQsolution2_+YQsolution1_)/2.0;
            rim_radius_X_ = sqrt(asquared_*(1.0-((((YQsolution1_+YQsolution2_)/2.0)-y0_)^2.0)/bsquared_));
            rim_radius_Y_ = ((YQsolution1_-YQsolution2_)/2.0);
            % Offset Paraboloid PROPERTY ASSIGNMENT AREA END..
            
            % Offset Paraboloid OBJECT CONSTUCTION AREA START.
            MRParaboloid = paraboloid([0;0;0],obj.F);
            MRRim = ellipticalRim([rim_move_X_;0],[rim_radius_X_;rim_radius_Y_],0);
            MRCoor = coordinateSystem([0,0,0]);
            
            obj.MR = reflector(MRParaboloid,MRRim,MRCoor);
            
            obj.coor_Feed = coordinateSystem([0,0,obj.F],[1;0;0],[0;1;0],MRCoor); %
            obj.coor_Feed = obj.coor_Feed.rotGRASP([((180.0-obj.beta/pi*180.0)*pi/180),0,0]);
            % Offset Paraboloid OBJECT CONSTUCTION AREA END.
        end
        
        function obj = translateReflectorCoordinateSystem(obj,delta)
            obj.MR.coor = obj.MR.coor.translate(delta);
        end
        
        function obj = rotGRASPReflectorCoordinateSystem(obj,angGRASP)
            obj.MR.coor = obj.MR.coor.rotGRASP(angGRASP);
        end
        
        function rho = getThRhMapping(obj,th)
            
        end
        
        function pathLengthStruct = getPathLength(obj,th,ph)
            % Calculate path length from:
            % feed to PR: FP
            % from PR to aperture: PA
            % from feed to aperture: FA
        end
        
        %% Plotting
        %PLOT NOT IMPLEMENTED.
        function plot(obj,N)
            if nargin == 1
                N = 101;
            end
            x = linspace(-obj.D/2,obj.D/2,N);
            y = zeros(size(x));
            
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
            
            obj.MR.plot(N);
            
            %obj.feedCoor.plot(obj.D/10);
            
            if coorFlag(1), obj.coor_Feed.plot(); end
            %if coorFlag(2), obj.apertureCoor.plot(obj.D/10); end
            if coorFlag(2), globalCoordSys = coordinateSystem(); globalCoordSys.plot(); end
            if coorFlag(3), obj.MR.coor.plot(); end
            xlabel('x-axis (m)')
            ylabel('y-axis (m)')
            zlabel('z-axis (m)')
        end
        
    end
    
end