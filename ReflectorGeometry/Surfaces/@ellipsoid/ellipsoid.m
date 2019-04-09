classdef ellipsoid
    % Rotationally symmetric ellipsoid (spheroid) - two axis have the same
    % length. See the 2002 Granet paper and GRASP technical description for details. 
    
    properties (SetAccess = private)
        vertexDistance = 2 % Major axis length in (m) = 2a
        fociDistance = 1 % Interfocal distance in (m) = 2f
        a % Vertex (Major) half distance
        f % Focus half distance - called c in the GRASP technical description and f in the Granet paper
        e % Eccentricity
        b % Minor half axis
    end
    
    methods
        function obj = ellipsoid(vertexDistance,fociDistance)
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
            obj.b = sqrt(obj.a.^2 - obj.f.^2);
        end
        
        function z = getZ(obj,x,y)
            % Only return the right half of the surface - see Granet and
            % GRASP technical manual
            z = obj.a.*sqrt(1 - obj.getRho(x,y).^2./obj.b.^2) - obj.f;
        end
        
        function rho = getRho(obj,x,y)
            rho = hypot(x,y);
            % Limit rho to always provide a structure that does not close
            % on itself...
            rho(rho > obj.b) = obj.b; 
        end
        
        function r = getR(obj,x,y)
            % See the GRASP technical description CH6 for details
            % The right hand focus is at the origin here, so r and v is
            % associated with the origin focus, and rp and vp with the
            % second focus.  v and vp are in radians
            v = obj.getV(x,y);
            r = obj.a.*(1 - obj.e.^2)./(1 + obj.e.*cos(v));
        end
        
        function rp = getRp(obj,x,y)
            % See the GRASP technical description CH6 for details
            % The right hand focus is at the origin here, so r and v is
            % associated with the origin focus, and rp and vp with the
            % second focus.  v and vp are in radians
            vp = obj.getVp(x,y);
            rp = obj.a.*(1 - obj.e.^2)./(1 - obj.e.*cos(vp));
        end
        
        function v = getV(obj,x,y)
            % See the GRASP technical description CH6 for details
            % The right hand focus is at the origin here, so r and v is
            % associated with the origin focus, and rp and vp with the
            % second focus.  v and vp are in radians
            v = atan(obj.getRho(x,y)./obj.getZ(x,y));
        end

        function vp = getVp(obj,x,y)
            % See the GRASP technical description CH6 for details
            % The right hand focus is at the origin here, so r and v is
            % associated with the origin focus, and rp and vp with the
            % second focus.  v and vp are in radians
            v = obj.getV(x,y);
            vp = 2*atan((1 - obj.e)./(1 + obj.e).*tan(v./2));
        end
        
        function u = getU(obj,x,y)
            % See the GRASP technical description CH6 for details
            % The right hand focus is at the origin here, so r and v is
            % associated with the origin focus, and rp and vp with the
            % second focus.  v and vp are in radians
            v = obj.getV(x,y);
            u = (v - obj.getVp(v))./2;
        end
        
        function n = getNorm(obj,x,y)
            nx = 2.*x./obj.b.^2;
            ny = 2.*y./obj.b.^2;
            nz = 2.*(obj.getZ(x,y) + obj.f)./obj.a.^2;
            n = [nx(:),ny(:),nz(:)].';
            nMag = sqrt(sum(n.^2));
            n = -bsxfun(@rdivide,n,nMag);
        end
        
        function [Cz,Ct] = getCurvature(obj,x,y)
            % Cz is the curvature in the plane containing the z-axis, and Ct
            % is in the orthogonal plane
            v = obj.getV(x,y);
            Cz = obj.a./obj.b.^2.*cos(obj.getU(v)).^3;
            Ct = obj.a./obj.b.^2.*cos(obj.getU(v));
        end
        
    end
    
end