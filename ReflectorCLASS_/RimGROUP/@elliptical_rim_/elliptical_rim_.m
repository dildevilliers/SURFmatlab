classdef elliptical_rim_
   properties
      Value
   end
   methods
       function obj = elliptical_rim_(val)
           if nargin > 0
               if isnumeric(val)
                   obj.Value = val;
               else
                   error('Value must be numeric')
               end
           end
       end
       function r = roundOff(obj)
           r = round([obj.Value],2);
       end
       function r = multiplyBy(obj,n)
           r = [obj.Value] * n;
       end
       
   end
end