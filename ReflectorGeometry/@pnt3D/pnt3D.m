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
            % The arrays X, Y, and Z must be the same size (or any of them can be scalar).
            if nargin == 3
                % Get all the same size
                obj.x = (Y+eps)./(Y+eps).*(Z+eps)./(Z+eps).*X;
                obj.y = (X+eps)./(X+eps).*(Z+eps)./(Z+eps).*Y;
                obj.z = (X+eps)./(X+eps).*(Y+eps)./(Y+eps).*Z;
            end
            [obj.ph,obj.el,obj.r] = cart2sph(obj.x,obj.y,obj.z);
            obj.th = pi/2 - obj.el;
            obj.rho = hypot(obj.x,obj.y);
        end
        
        function S = size(obj)
            S = size(obj.x);
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
            
            parse(parseobj, obj, varargin{:});

            marker = parseobj.Results.marker;
            markerEdgeColor = parseobj.Results.markerEdgeColor;
            markerFaceColor = parseobj.Results.markerFaceColor;
            markerSize = parseobj.Results.markerSize;

            plot3(obj.x(:),obj.y(:),obj.z(:),'marker',marker,'markerEdgeColor',markerEdgeColor,'markerFaceColor',markerFaceColor,'markerSize',markerSize), grid on
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