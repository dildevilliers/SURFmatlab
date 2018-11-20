function obj = setYrange(obj,type)
% Attempts to set the y-range (th, el, or al) for the angular
% gridTypes in the FarField object
% type = 180:
% [0,180] gridType = 'PhTh'
% [-90,90] gridType = 'AzEl' | 'ElAz'
% type = 360:
% [0,360] for xRangeType = 'pos'
% [-180,180] for xRangeType = 'sym'
% Redundant fields are replaced by valid ones as far as possible
% The resulting object has the same number
% of field points as the input object.

assert(type == 180 || type == 360,'Unknown type: Should be 180 or 360');

if strcmp(obj.gridType,'PhTh')
    if type == 360
        if strcmp(obj.xRangeType,'pos')
            % ToDo
        elseif strcmp(obj.xRangeType,'sym')
            % ToDo
        end
    elseif type == 180
        % ToDo
    end
elseif strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
    % ToDo
else
    warning(['Cant shift a polar grid like ', obj.gridType, ' on a cartesian grid']);
end

end