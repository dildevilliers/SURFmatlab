classdef paraboloid
   properties
      vertex(1,1) pnt3D = pnt3D(0,0,0)
      focalLength(1,1) double {mustBeReal, mustBeFinite} = 1 % in (m) 
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
       end
       
       function z = getZ(obj,x,y)
           z = ((x - obj.vertex.x).^2 + (y - obj.vertex.y).^2)./(4.*obj.focalLength) + obj.vertex.z;
       end
       
       function n = getNorm(obj,x,y)
           nx = (x - obj.vertex.x)./(2*obj.focalLength);
           ny = (y - obj.vertex.y)./(2*obj.focalLength);
           nz = ones(size(nx)).*(-1);  
           n = [nx(:),ny(:),nz(:)].';
           nMag = sqrt(sum(n.^2));
           n = -bsxfun(@rdivide,n,nMag);
       end
       
       function C = getCurvature(obj,x,y)
           
       end
           
   end
end