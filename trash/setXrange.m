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

mustBeMember(type,{'pos','sym'})
if strcmp(type,'pos')
    t = 'p';
    if strcmp(obj.xRangeType,'pos'), return; end
elseif strcmp(type,'sym')
    t = 's';
    if strcmp(obj.xRangeType,'sym'), return; end
end
% tol = 1e-10;
tol = 10^(-obj.nSigDig);
if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
    if t == 'p'
%         iout = find(obj.x == -pi);   % Redundant
%         iin = find(obj.x == 0);      % Will become redundant after inserting
        iout = find(abs(obj.x + pi) < tol);   % Redundant
        iin = find(abs(obj.x - 0) < tol);      % Will become redundant after inserting
        if strcmp(obj.xRangeType,'pos'), return; end
    elseif t == 's'
%         iout = find(obj.x == 2*pi);   % Redundant
%         iin = find(obj.x == pi);      % Will become redundant after inserting
        iout = find(abs(obj.x - 2*pi) < tol);   % Redundant
        iin = find(abs(obj.x - pi) < tol);      % Will become redundant after inserting
        if strcmp(obj.xRangeType,'sym'), return; end
    end
    Nredun = min(numel(iout),numel(iin));    % How many we have to remove/add
    % Truncate the indexes to the shortest length
    iout = iout(1:Nredun);
    iin = iin(1:Nredun);
    yAdd = obj.y(iin);
    E1Add = obj.E1(iin,:);
    E2Add = obj.E2(iin,:);
    redunFound = Nredun > 0;
    % First remove the ph=-180 from the matrix...
    if redunFound
        obj.x(iout) = [];
        obj.y(iout) = [];
        obj.E1(iout,:) = [];
        obj.E2(iout,:) = [];
    end
    % Now shift the x values
    if t == 'p'
        outOfRangeInd = find(obj.x < 0);
        obj.x(outOfRangeInd) = obj.x(outOfRangeInd) + 2*pi;
    elseif t == 's'
        outOfRangeInd = find(obj.x > pi);
        obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - 2*pi;
    end
    % if the -180 was removed, add a 360 cut
    if redunFound
        if t == 'p'
            xAdd = ones(numel(iin),1).*2.*pi;
        elseif t == 's'
            xAdd = ones(numel(iin),1).*-pi;
        end
        obj.x = [obj.x;xAdd(1:Nredun)];
        obj.y = [obj.y;yAdd(1:Nredun)];
        obj.E1 = [obj.E1;E1Add(1:Nredun,:)];
        obj.E2 = [obj.E2;E2Add(1:Nredun,:)];
    end
    % Sort
    obj = obj.sortGrid;
else
    warning(['Cant shift a polar grid like ', obj.gridType, ' on a cartesian grid']);
end

end