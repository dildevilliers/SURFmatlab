classdef ellipsoid
    % Rotationally symmetric ellipsoid (spheroid) - two axis have the same
    % length. See the 2002 Granet paper and GRASP technical description for details. 
    
    properties (SetAccess = private)
        vertexDistance = 2 % Major axis length in (m) = 2a
        fociDistance = 1 % Interfocal distance in (m) = 2f
        rotAng = 0     % Rotation angle in rad from negative z-axis
        a % Vertex (Major) half distance
        f % Focus half distance - called c in the GRASP technical description and f in the Granet paper
        e % Eccentricity
        b % Minor half axis
        F1 % First focus position (origin of coor)
        F0 % Second focus position (where the feed goes)
    end
    
    properties (SetAccess = private, Hidden = true)
        % Az^2 + Bxz + Cz + Dx^2 + Ex + Fy^2 - G = 0
        A
        B
        C
        D
        E
        F
        G
        coorRot
    end
    
    methods
        function obj = ellipsoid(vertexDistance,fociDistance,rotAng)
            if nargin >= 2
                obj.vertexDistance = vertexDistance;
                obj.fociDistance = fociDistance;
            end
            if nargin >= 3, obj.rotAng = rotAng; end
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
            % Focus positions
            coorBase = CoordinateSystem();
            obj.coorRot = coorBase.rotY(-obj.rotAng);
            obj.F1 = Pnt3D(0,0,0);
            obj.F0 = Pnt3D(0,0,2*obj.f);
            obj.F0 = obj.F0.changeBase(coorBase,obj.coorRot);
            % Plenty of algebra to get these...
            th = obj.rotAng;
            obj.A = sin(th).^2./obj.b.^2 + cos(th).^2./obj.a.^2;
            obj.D = cos(th).^2./obj.b.^2 + sin(th).^2./obj.a.^2;
            obj.B = (obj.e./obj.b).^2.*sin(2.*th);
            obj.C = -2.*(obj.e./obj.a).*cos(th);
            obj.E = +2.*(obj.e./obj.a).*sin(th);
            obj.F = 1./obj.b.^2;
            obj.G = (obj.b./obj.a).^2;
        end
        
        function [zn,zp] = getZ(obj,x,y)
            % Only return the negative z of the surface 
            % see GRASP technical manual Figure 2-3
            Alp = getAlp(obj);
            Bet = getBet(obj,x);
            Gam = getGam(obj,x,y);
            Det = Bet.^2 - 4.*Alp.*Gam;
            zn = -(Bet + sqrt(Det))./(2.*Alp);
            zp = -(Bet - sqrt(Det))./(2.*Alp);
            zp(Det < 0) = max(max(zp(Det > 0)));  % Do something about the points outside the acceptable range
            zn(Det < 0) = max(max(zn(Det > 0)));  % Do something about the points outside the acceptable range
        end
        
        function xz = getXZ(obj,N)
            % Returns a vector of length N of all the points around
            % the ellipse in the x-z plane.  Handy for plotting
            if nargin < 2, N = 200; end
            % First get the points of a non-rotated version
            ell0 = ellipsoid(obj.vertexDistance,obj.fociDistance,0);
            xmin0 = -ell0.b*0.999;
            xmax0 = ell0.b*0.999;
            x0 = linspace(xmin0,xmax0,round(N/2));
            y0 = zeros(size(x0));
            [zn0,zp0] = ell0.getZ(x0,y0);
            x0vect = [x0,x0].';
            y0vect = [y0,y0].';
            z0vect = [zn0,zp0].';
            p0 = Pnt3D(x0vect,y0vect,z0vect);
            xz = p0.changeBase(CoordinateSystem,obj.coorRot);
        end
        
        function alpha = getAlp(obj)
            alpha = obj.A;
        end
        
        function beta = getBet(obj,x)
            beta = obj.B.*x + obj.C;
        end
        
        function gamma = getGam(obj,x,y)
           gamma = obj.D.*x.^2 + obj.E.*x + obj.F.*y.^2 - obj.G; 
        end
        
        function rho = getRho(obj,x,y)
            % Points on surface in global coor
            surfPnt = Pnt3D(x,y,obj.getZ(x,y)); 
            % Get in the local rotated coordinate system
            surfPntLoc = surfPnt.changeBase(obj.coorRot,CoordinateSystem);
            rho = surfPntLoc.rho;
        end
        
        function r = getR(obj,x,y)
            % See the GRASP technical description CH6 for details
            surfPnt = Pnt3D(x,y,obj.getZ(x,y)); 
            rVect = surfPnt - obj.F0;
            r = rVect.r;
        end
        
        function rp = getRp(obj,x,y)
            % See the GRASP technical description CH6 for details
            surfPnt = Pnt3D(x,y,obj.getZ(x,y)); 
            rpVect = surfPnt - obj.F1;
            rp = rpVect.r;
        end
        
        function v = getV(obj,x,y)
            % See the GRASP technical description CH6 for details
            % Always returns the acute angle (not really used)
            % v and vp are in radians
            v = asin(obj.getRho(x,y)./obj.getR(x,y));
        end

        function vp = getVp(obj,x,y)
            % See the GRASP technical description CH6 for details
            % Always returns the acute angle (not really used)
            % v and vp are in radians
            vp = asin(obj.getRho(x,y)./obj.getRp(x,y));
        end
        
        function u = getU(obj,x,y)
            % See the GRASP technical description CH6 for details
            % Use to cosine rule to get it directly from r and rp
            r = getR(obj,x,y);
            rp = getRp(obj,x,y);
            u = 0.5.*acos((r.^2 + rp.^2 - (2.*obj.f).^2)./(2.*r.*rp));
        end
        
        function n = getNorm(obj,x,y)
            % Get the analytical derivative from the hand calculated
            % formula
            z = obj.getZ(x,y);
            nx = obj.B.*z + 2.*obj.D.*x + obj.E;
            ny = 2.*obj.F.*y;
            nz = 2.*obj.A.*z + obj.B.*x + obj.C;
%             nx = 2.*x./obj.b.^2;
%             ny = 2.*y./obj.b.^2;
%             nz = 2.*(obj.getZ(x,y) + obj.f)./obj.a.^2;
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