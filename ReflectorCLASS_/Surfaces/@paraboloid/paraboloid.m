classdef paraboloid
   properties
      vertex(3,1) double {mustBeReal, mustBeFinite} = [0;0;0] % [x,y,z] in (m) of vertex
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
           z = ((x - obj.vertex(1)).^2 + (y - obj.vertex(2)).^2)./(4.*obj.focalLength) + obj.vertex(3);
       end
   end
end