classdef ArrayElements
    %ARRAYELEMENTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        arrayPositions(:,1) pnt3D % {mustBeFinite}
        elementPatterns(:,1) FarField 
    end
    
    methods
        function obj = ArrayElements(arrayPositions,elementPatterns)
            obj.arrayPositions = arrayPositions;
            
            if nargin < 2
                FF = FarField;
                FF = FF.makeIsotropic();
                obj.elementPatterns(size(arrayPositions),1) = FF;
            else
                if size(arrayPositions) ~= size(elementPatterns)
                   error('The elementPatterns must be equal in size to the arrayPositions');
                end
                obj.elementPatterns = elementPatterns;
            end

        end
    end
end

