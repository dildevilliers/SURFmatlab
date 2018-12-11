function obj = setYrange(obj,type)
% Attempts to set the y-range (th, el, or al) for the angular
% gridTypes in the FarField object
% type = 180:
% [0,180] gridType = 'PhTh'
% [-90,90] gridType = 'AzEl' | 'ElAz'
% type = 360:
% [0,360] for xRangeType = 'pos' - xRange = [0,180]
% [-180,180] for xRangeType = 'sym' - xRange = [-90,90]
% Redundant fields are replaced by valid ones as far as possible
% The resulting object has the same number
% of field points as the input object.

assert(type == 180 || type == 360,'Unknown type: Should be 180 or 360');

iout = [];
iin = [];

if strcmp(obj.gridType,'PhTh')
    if type == 360
        if strcmp(obj.xRangeType,'pos')
            % Redundant pole: remove ph > 180, th = 180; add ph = 0, th = 180:360
            iout = find(obj.x  > pi & abs(obj.y - pi) < 1e-10);   % Redundant
            iin = find(abs(obj.x - 0) < 1e-10);                   % Will become redundant after inserting
            [obj,iin] = removeRedun(obj,iout,iin);
            outOfRangeInd = find(obj.x > pi);
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - pi;
            obj.y(outOfRangeInd) = 2*pi - obj.y(outOfRangeInd);
            % Insert new redundant points
            obj.x = [obj.x;ones(numel(iin),1).*0];
            obj.y = [obj.y;2*pi - obj.y(iin)];
        elseif strcmp(obj.xRangeType,'sym')
%             Redundant pole: remove |ph| > 90, th = 0 and  ph = -180, th = 0:180; add ph = abs(90), th = 0:180
            iout = unique([find((obj.x > pi/2 | obj.x < -pi/2) & abs(obj.y - 0) < 1e-10); find(abs(obj.x + pi) < 1e-10)]);   % Redundant
            iin1 = find(abs(obj.x + pi/2) < 1e-10);
            iin2 = find(abs(obj.x - pi/2) < 1e-10);
            iin = unique([iin1;iin2]);
            [obj,iin] = removeRedun(obj,iout,iin);
            outOfRangeInd = find(obj.x > pi/2 | obj.x < -pi/2);
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - sign(obj.x(outOfRangeInd)).*pi;
            obj.y(outOfRangeInd) = -obj.y(outOfRangeInd);
            % Insert (as many as possible) new redundant points
            niinKeep = numel(iin) - numel(iin1);
            xAdd = [ones(numel(iin1),1).*-pi/2;ones(niinKeep,1).*pi/2];
            yAdd = [obj.y(iin1);-obj.y(iin2(1:niinKeep))];
            obj.x = [obj.x;xAdd];
            obj.y = [obj.y;yAdd];
        end
    elseif type == 180
        % ToDO: Redundancies...
        if strcmp(obj.xRangeType,'pos')
            outOfRangeInd = find(obj.y > pi);
            obj.y(outOfRangeInd) = 2*pi - obj.y(outOfRangeInd);
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) + pi;
        elseif strcmp(obj.xRangeType,'sym')
            outOfRangeInd = find(obj.y < 0);
            obj.y(outOfRangeInd) = -obj.y(outOfRangeInd);
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - sign(obj.x(outOfRangeInd)).*pi;
        end
    end
elseif strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
    % ToDo
else
    warning(['Cant shift a polar grid like ', obj.gridType, ' on a cartesian grid']);
end
% Insert new redundant E-field values
obj.E1 = [obj.E1;obj.E1(iin,:)];
obj.E2 = [obj.E2;obj.E2(iin,:)];
obj.E3 = [obj.E3;obj.E3(iin,:)];
% Sort
obj = obj.sortGrid;
% Update the object descriptive parameters
obj = setRangeTypes(obj);
end


% Function to remove the redundant values, and return the redacted vector
% of insertion indexes
function [objSmall,iinSmall] = removeRedun(objFull,iout,iin)
objSmall = objFull;
Nredun = min(numel(iout),numel(iin));    % How many we have to remove/add
% Truncate the indexes to the shortest length
iout = iout(1:Nredun);
iinSmall = iin(1:Nredun);
redunFound = Nredun > 0;
% First remove the ph=-180 from the matrix...
if redunFound
    objSmall.x(iout) = [];
    objSmall.y(iout) = [];
    objSmall.E1(iout,:) = [];
    objSmall.E2(iout,:) = [];
    objSmall.E3(iout,:) = [];
end
end