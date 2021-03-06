classdef paraboloid
    % See the GRASP technical description Chapter 6 for details
   properties (SetAccess = private)
      vertex = Pnt3D(0,0,0)
      focalLength = 1 % in (m) 
      F0 = Pnt3D(0,0,1) % Focus point
   end
   methods
       function obj = paraboloid(ver,F)
           if nargin == 0
           elseif nargin == 1
               obj.vertex = ver;
           else
               obj.vertex = ver;
               obj.focalLength = F;
           end
           obj = obj.setParams;
       end
       
       function obj = setParams(obj)
            obj.F0 = obj.vertex + Pnt3D(0,0,obj.focalLength);
       end
       
       function obj = setVertex(obj,vert)
           obj.vertex = vert;
           obj = obj.setParams;
       end
       
       function obj = setFocalLength(obj,focLen)
           obj.focalLength = focLen;
           obj = obj.setParams;
       end
       
       function rho = getRho(obj,x,y)
           rho = hypot(x-obj.vertex.x,y-obj.vertex.y);
       end
       
       function z = getZ(obj,x,y)
           z = (obj.getRho(x,y).^2)./(4.*obj.focalLength) + obj.vertex.z;
       end
       
       function r = getR(obj,x,y)
           r = obj.focalLength./cos(obj.getU(x,y)).^2;
       end
       
       function u = getU(obj,x,y)
           rho = obj.getRho(x,y);
           u = atan(rho./(2*obj.focalLength));
       end
       
       function n = getNorm(obj,x,y)
           nx = (x - obj.vertex.x)./(2*obj.focalLength);
           ny = (y - obj.vertex.y)./(2*obj.focalLength);
           nz = ones(size(nx)).*(-1);  
           n = [nx(:),ny(:),nz(:)].';
           nMag = sqrt(sum(n.^2));
           n = -bsxfun(@rdivide,n,nMag);
       end
       
       function [Cz,Ct] = getCurvature(obj,x,y)
           % Cz is the curvature in the plane containing the z-axis, and Ct
           % is in the orthogonal plane
           Cz = cos(obj.getU(x,y)).^3./(2.*obj.focalLength);
           Ct = cos(obj.getU(x,y))./(2.*obj.focalLength);
       end
           
   end
end