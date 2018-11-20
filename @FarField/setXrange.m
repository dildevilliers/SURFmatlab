function obj = setXrange(obj,type)
% Attempts to set the x-range (ph, az, or ep) for the angular
% gridTypes in the FarField object
% type = 'pos':
% [-180,180] is transformed to [0,360]
% with the redundant -180 cut replaced by a redundant 360 cut.
% type = 'sym':
% [0,360] is transformed to [-180,180]
% with the redundant 360 cut replaced by a redundant -180 cut.
% The resulting object has the same number
% of field points as the input object.

assert((strcmp(type,'pos') || strcmp(type,'sym')),'Unknown type: Should be pos or sym');
if strcmp(type,'pos')
    t = 'p';
elseif strcmp(type,'sym')
    t = 's';
end

if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
    if t == 'p'
        iout = find(obj.x == -pi);   % Redundant
        iin = find(obj.x == 0);      % Will become redundant
    elseif t == 's'
        iout = find(obj.x == 2*pi);   % Redundant
        iin = find(obj.x == pi);      % Will become redundant
    end
    Nredun = min(numel(iout),numel(iin));    % How many we have to remove/add
    % Truncate the indexes to the shortest length
    iout = iout(1:Nredun);
    iin = iin(1:Nredun);
    redunFound = Nredun > 0;
    % First remove the ph=-180 from the matrix...
    if redunFound
        obj.x(iout) = [];
        obj.y(iout) = [];
        obj.E1(iout,:) = [];
        obj.E2(iout,:) = [];
        obj.E3(iout,:) = [];
    end
    % Now shift the x values
    if t == 'p'
        rangeEdge = find(obj.x < 0);
        obj.x(rangeEdge) = obj.x(rangeEdge) + 2*pi;
    elseif t == 's'
        rangeEdge = find(obj.x > pi);
        obj.x(rangeEdge) = obj.x(rangeEdge) - 2*pi;
    end
    % if the -180 was removed, add a 360 cut
    if redunFound
        if t == 'p'
            obj.x = [obj.x;ones(numel(iin),1).*2.*pi];
        elseif t == 's'
            obj.x = [obj.x;ones(numel(iin),1).*-pi];
        end
        obj.y = [obj.y;obj.y(iin)];
        obj.E1 = [obj.E1;obj.E1(iin,:)];
        obj.E2 = [obj.E2;obj.E2(iin,:)];
        obj.E3 = [obj.E3;obj.E3(iin,:)];
    end
    % Sort
    obj = obj.sortGrid;
else
    warning(['Cant shift a polar grid like ', obj.gridType, ' on a cartesian grid']);
end
obj = setRangeTypes(obj);
end