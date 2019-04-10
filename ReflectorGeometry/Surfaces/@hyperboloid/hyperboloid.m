classdef hyperboloid
    % Rotationally symmetric hyperboloid (spheroid).
    % See the 2002 Granet paper and GRASP technical description for details. 
    
    properties (SetAccess = private)
        vertexDistance = 1 % Vertex seperation distance in (m) = 2a - positive is convex (like cassegrain), negative is concave (compact range)
        fociDistance = 2 % Interfocal distance in (m) = 2f
        a % Vertex seperation half distance and sign
        f % Focus half distance - called c in the GRASP technical description and f in the Granet paper
        e % Eccentricity
        b 
        F1 % First focus position (origin of coor)
        F0 % Second focus position (where the feed goes)
    end
    
    methods
        function obj = hyperboloid(vertexDistance,fociDistance)
            if nargin == 0
            elseif nargin == 2
                obj.vertexDistance = vertexDistance;
                obj.fociDistance = fociDistance;
            end
            obj = obj.setParams;
        end
        
        function obj = setVertexDistance(obj,vD)
            obj.vertexDistance = vD;
            obj = obj.setParams;
        end
        
        function obj = setFociDistance(obj,fD)
            obj.vertexDistance = fD;
            obj = obj.setParams;
        end
        
        function obj = setParams(obj)
            obj.a = obj.vertexDistance/2;
            obj.f = obj.fociDistance/2;
            obj.e = obj.f./obj.a;
            obj.b = sqrt(obj.f.^2 - obj.a.^2);
            obj.F1 = pnt3D(0,0,0);
            obj.F0 = pnt3D(0,0,-2*obj.f);
        end
        
        function z = getZ(obj,x,y)
            % Only return the right half of the surface - see Granet and
            % GRASP technical manual
            z = obj.a.*sqrt(1 + obj.getRho(x,y).^2./obj.b.^2) - obj.f;
        end
        
        function rho = getRho(obj,x,y)
            rho = hypot(x,y);
        end
        
        function r = getR(obj,x,y)
            % See the GRASP technical description CH6 for details
            % The right hand focus is at the origin here, so r and v is
            % associated with the origin focus, and rp and vp with the
            % second focus.  v and vp are in radians
            v = obj.getV(x,y);
            r = obj.a.*(obj.e.^2 - 1)./(obj.e.*cos(v) + 1);
        end
        
        function rp = getRp(obj,x,y)
            % See the GRASP technical description CH6 for details
            % The right hand focus is at the origin here, so r and v is
            % associated with the origin focus, and rp and vp with the
            % second focus.  v and vp are in radians
            vp = obj.getVp(x,y);
            rp = obj.a.*(obj.e.^2 - 1)./(obj.e.*cos(vp) - 1);
        end
        
        function v = getV(obj,x,y)
            % See the GRASP technical description CH6 for details
            % The right hand focus is at the origin here, so r and v is
            % associated with the origin focus, and rp and vp with the
            % second focus.  v and vp are in radians
            v = -atan(obj.getRho(x,y)./obj.getZ(x,y));
        end

        function vp = getVp(obj,x,y)
            % See the GRASP technical description CH6 for details
            % The right hand focus is at the origin here, so r and v is
            % associated with the origin focus, and rp and vp with the
            % second focus.  v and vp are in radians
            v = obj.getV(x,y);
            vp = 2*atan((obj.e - 1)./(obj.e + 1).*tan(v./2));
        end
        
        function u = getU(obj,x,y)
            % See the GRASP technical description CH6 for details
            % The right hand focus is at the origin here, so r and v is
            % associated with the origin focus, and rp and vp with the
            % second focus.  v and vp are in radians
            v = obj.getV(x,y);
            u = (v + obj.getVp(v))./2;
        end
        
        function n = getNorm(obj,x,y)
            nx = -2.*x./obj.b.^2;
            ny = -2.*y./obj.b.^2;
            nz = 2.*(obj.getZ(x,y) + obj.f)./obj.a.^2;
            n = [nx(:),ny(:),nz(:)].';
            nMag = sqrt(sum(n.^2));
            n = -bsxfun(@rdivide,n,nMag);
        end
        
        function [Cz,Ct] = getCurvature(obj,x,y)
            % Cz is the curvature in the plane containing the z-axis, and Ct
            % is in the orthogonal plane
            Cz = obj.a./obj.b.^2.*cos(obj.getU(x,y)).^3;
            Ct = obj.a./obj.b.^2.*cos(obj.getU(x,y));
        end
        
    end
    
end