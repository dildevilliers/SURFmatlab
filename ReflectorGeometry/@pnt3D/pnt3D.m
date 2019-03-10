classdef pnt3D
    properties
        x double {mustBeReal, mustBeFinite} = 1 % in (m)
        y double {mustBeReal, mustBeFinite} = 1 % in (m)
        z double {mustBeReal, mustBeFinite} = 1 % in (m)
    end
    
    properties (SetAccess = private)
        % Values in other coordinate systems
        th % polar angle in radians
        ph % azimuth angle in radians
        el % Elevation angle in radians
        r  % distance from origin
        rho % distance from z-axis
    end
    
    methods
        function obj = pnt3D(X,Y,Z)
            if nargin == 3
                % Get all the same size
                obj.x = (Y+eps)./(Y+eps).*(Z+eps)./(Z+eps).*X;
                obj.y = (X+eps)./(X+eps).*(Z+eps)./(Z+eps).*Y;
                obj.z = (X+eps)./(X+eps).*(Y+eps)./(Y+eps).*Z;
            end
            obj = obj.setProps;
        end
        
        function obj = setProps(obj)
            [obj.ph,obj.el,obj.r] = cart2sph(obj.x,obj.y,obj.z);
            obj.th = pi/2 - obj.el;
            obj.rho = hypot(obj.x,obj.y);
        end
        
        function S = size(obj)
            S = size(obj.x);
        end
        
        function obj = getNpts(obj,I)
            % Returns only the points in the indexes I in the new object
            obj = pnt3D(obj.x(I),obj.y(I),obj.z(I));
        end
        
        function B = isequal(obj1,obj2)
            tol = eps;
            D = obj1-obj2;
            B = all(all([abs(D.x),abs(D.y),abs(D.z)] < tol));
        end
        
        function obj = plus(obj,obj2)
            obj.x = obj.x+obj2.x;
            obj.y = obj.y+obj2.y;
            obj.z = obj.z+obj2.z;
            obj = obj.setProps;
        end
        
        function obj = minus(obj,obj2)
            obj.x = obj.x-obj2.x;
            obj.y = obj.y-obj2.y;
            obj.z = obj.z-obj2.z;
            obj = obj.setProps;
        end
        
        function D = distanceCart(obj,obj2)
            objD = obj - obj2;
            D = objD.r;
        end
        
        function obj = translate(obj,DEL)
            % translates the point(s) by DEL = [dx;dy;dz]
            Xp = obj.pointMatrix + DEL;
%             Xp = trans3D(pointMatrix(obj),DEL(:));
            obj.x = reshape(Xp(1,:),size(obj));
            obj.y = reshape(Xp(2,:),size(obj));
            obj.z = reshape(Xp(3,:),size(obj));
            obj = obj.setProps;
        end
        
        function X = pointMatrix(obj)
            % Rows ar [x;y;z]
            X = [obj.x(:),obj.y(:),obj.z(:)].';
        end
        
        
        function obj = changeBase(obj,coor_new,coor_base)
            % Transforms the points from the coordinate system coor_base,
            % to the new coordinate system coor_new through translation and rotation.
            
            if nargin == 2
                coor_base = coordinateSystem();
            end
            % Move points to new coordinate origin reference 
            U = pointMatrix(obj) - coor_new.origin.pointMatrix;
            % Rotate the points in the origin reference
            Q = dirCosine(coor_new,coor_base);
            Uprime = Q\U;
            % Move to new coordinate base 
            Uprime = Uprime + coor_base.origin.pointMatrix;
            % Make the object
            obj.x = reshape(Uprime(1,:),size(obj));
            obj.y = reshape(Uprime(2,:),size(obj));
            obj.z = reshape(Uprime(3,:),size(obj));
            obj = obj.setProps;
        end
        
        function plot(obj,varargin)
            parseobj = inputParser;
            parseobj.FunctionName = 'plot';
            
            typeValidationObj = @(x) validateattributes(x,{'pnt3D'},{'numel',1},'plot','obj',1);
            addRequired(parseobj,'obj',typeValidationObj);
            
            typeValidationMarker = @(x) validateattributes(x,{'char'},{},'plot','marker');
            addParameter(parseobj,'marker','.',typeValidationMarker);
            
            typeValidationMarkerEdgeColor = @(x) validateattributes(x,{'char','double'},{},'plot','markerEdgeColor');
            addParameter(parseobj,'markerEdgeColor','k',typeValidationMarkerEdgeColor);
            
            typeValidationMarkerFaceColor = @(x) validateattributes(x,{'char','double'},{},'plot','markerFaceColor');
            addParameter(parseobj,'markerFaceColor','none',typeValidationMarkerFaceColor);
            
            typeValidationMarkerSize = @(x) validateattributes(x,{'double'},{'real','positive'},'plot','markerSize');
            addParameter(parseobj,'markerSize',10,typeValidationMarkerSize);
            
            typeValidationLineStyle = @(x) validateattributes(x,{'char'},{},'plot','lineStyle');
            addParameter(parseobj,'lineStyle','none',typeValidationLineStyle);
            
            typeValidationLineColor = @(x) validateattributes(x,{'char','double'},{},'plot','lineColor');
            addParameter(parseobj,'lineColor','k',typeValidationLineColor);
            
            typeValidationLineWidth = @(x) validateattributes(x,{'double'},{'real','positive'},'plot','lineWidth');
            addParameter(parseobj,'lineWidth',1,typeValidationLineWidth);
            
            parse(parseobj, obj, varargin{:});
            
            marker = parseobj.Results.marker;
            markerEdgeColor = parseobj.Results.markerEdgeColor;
            markerFaceColor = parseobj.Results.markerFaceColor;
            markerSize = parseobj.Results.markerSize;
            lineStyle = parseobj.Results.lineStyle;
            lineColor = parseobj.Results.lineColor;
            lineWidth = parseobj.Results.lineWidth;
            
            plot3(obj.x(:),obj.y(:),obj.z(:),'linestyle',lineStyle,...
                'color',lineColor,'lineWidth',lineWidth,...
                'marker',marker,'markerEdgeColor',markerEdgeColor,...
                'markerFaceColor', markerFaceColor,'markerSize',markerSize), grid on
            xlabel('x (m)')
            ylabel('y (m)')
            zlabel('z (m)')
        end
    end
    
    methods (Static = true)
        function obj = sph(PH,TH,R)
            % Define in spherical coordinates
            [X,Y,Z] = sph2cart(PH,pi/2 - TH,R);
            obj = pnt3D(X,Y,Z);
        end
        
        function obj = pol(PH,RHO,Z)
            % Define in polar coordinates
            if nargin == 2
                Z = 0;
            end
            [X,Y] = pol2cart(PH,RHO,Z);
            obj = pnt3D(X,Y,Z);
        end
    end
    
end