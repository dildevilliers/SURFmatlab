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

[iout,iin,xAdd,yAdd,E1Add,E2Add,E3Add] = deal([]);
Nredun = 0;

% eps = 1e-10;
eps = 1e-1*min([diff(unique(obj.x));diff(unique(obj.y))]);

% For all cases some redundant points will be removed, and some new ones
% inserted into the grid. It depends on how the space is cut and rotated
% where to put in and take out points...  Done now on a case-by-case basis
% because I can't think of a general algorithm yet
if strcmp(obj.gridType,'PhTh')
    if type == 360
        if strcmp(obj.xRangeType,'pos')
            iout = find(obj.x  > (pi+eps) & abs(obj.y - pi) < eps);   % Redundant
            iin = find(abs(obj.x - 0) < eps);                   % Will become redundant after inserting
            xAdd = obj.x(iin);
            yAdd = 2*pi - obj.y(iin);
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.x > (pi+eps));
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - pi;
            obj.y(outOfRangeInd) = 2*pi - obj.y(outOfRangeInd);
        elseif strcmp(obj.xRangeType,'sym')
            iout = [find((obj.x >= (-pi/2-eps) & obj.x <= (pi/2+eps)) & abs(obj.y - 0) < eps); find(abs(obj.x + pi) < eps)];   % Redundant
            iin = find(abs(obj.x + pi/2) < eps | abs(obj.x - pi/2) < eps);
            xAdd = obj.x(iin) - sign(obj.x(iin)).*pi;
            yAdd = -obj.y(iin);
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.x > (pi/2+eps) | obj.x < (-pi/2-eps));
            signPos = sign(obj.x(outOfRangeInd));
            signPos(signPos == 0) = 1;
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - signPos.*pi;
            obj.y(outOfRangeInd) = -obj.y(outOfRangeInd);
        end
    elseif type == 180
        if strcmp(obj.xRangeType,'pos')
            iout = find(obj.y <= (pi+eps) & abs(obj.x - pi) < eps);   % Redundant
            iin = find(obj.x <= (pi+eps) & abs(obj.y - pi) < eps);    % Will become redundant after inserting
            xAdd = obj.x(iin) + pi;
            yAdd = 2*pi - obj.y(iin);
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.y > (pi+eps));
            obj.y(outOfRangeInd) = 2*pi - obj.y(outOfRangeInd);
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) + pi;
        elseif strcmp(obj.xRangeType,'sym')
%             iout = unique(find((abs(obj.x - pi/2) < eps | abs(obj.x + pi/2) < eps) & (obj.y <= eps)));   % Redundant
            iout = find((abs(obj.x - pi/2) < eps | abs(obj.x + pi/2) < eps) & (obj.y <= eps));   % Redundant
            iin1 = find(abs(obj.y - 0) < eps);
            iin2 = find(abs(obj.x - 0) < eps & obj.y <= eps);
            iin = [iin1;iin2];
            signPosIn = sign(obj.x(iin1));
            signPosIn(signPosIn == 0) = 1;
            xAdd1 = obj.x(iin1) - signPosIn.*pi;
            xAdd2 = obj.x(iin2) + pi;
            xAdd = [xAdd1;xAdd2];
            yAdd = -obj.y(iin);
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.y < -eps);
            obj.y(outOfRangeInd) = -obj.y(outOfRangeInd);
            signPos = sign(obj.x(outOfRangeInd));
            signPos(signPos == 0) = 1;
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - signPos.*pi;
        end
    end
elseif strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
    if type == 360
        if strcmp(obj.xRangeType,'pos')
            iout = find(obj.x <= pi+eps & (abs(obj.y - pi/2) < eps | abs(obj.y + pi/2) < eps));   % Redundant
            iin1 = find(abs(obj.x - pi) < eps);                   % Will become redundant after inserting
            iin2 = find(abs(obj.y - 0) < eps);
            iin = [iin1;iin2];
            xAdd = [pi - obj.x(iin1);obj.x(iin2)];
            yAdd = [pi + obj.y(iin1);2*pi + obj.y(iin2)];
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift - first x > 180
            outOfRangeInd1 = find(obj.x > (pi+eps));
            obj.x(outOfRangeInd1) = obj.x(outOfRangeInd1) - pi;
            obj.y(outOfRangeInd1) = pi - obj.y(outOfRangeInd1);
            % Now for negative y
            outOfRangeInd2 = find(obj.y < -eps);
            obj.x(outOfRangeInd2) = obj.x(outOfRangeInd2);
            obj.y(outOfRangeInd2) = 2*pi + obj.y(outOfRangeInd2);
        elseif strcmp(obj.xRangeType,'sym')
            iout = find(abs(obj.x - pi) < eps | (obj.x <= (pi/2+eps) & obj.x >= (-pi/2-eps) & abs(obj.y - pi/2) < eps) | (obj.x <= (pi/2+eps) & obj.x >= (-pi/2-eps) & abs(obj.y + pi/2) < eps));   % Redundant
            iin1 = find(abs(obj.x + pi/2) < eps);
            iin2 = find(abs(obj.x - pi/2) < eps);
            iin3 = find(abs(obj.y - 0) < eps & (obj.x >= (pi/2-eps) | obj.x <= (-pi/2+eps)));
            xAdd1 = obj.x(iin1);
            xAdd2 = obj.x(iin2);
            xAdd3 = obj.x(iin3) - sign(obj.x(iin3)).*pi;    % Never zero so don't worry about making positive 
            xAdd = [xAdd1;xAdd2;xAdd3];
            iin12 = [iin1;iin2];
            signPos = sign(obj.y(iin12));
            signPos(signPos == 0) = 1;
            yAdd12 = [signPos.*pi - obj.y(iin12)];  % Include the -180 y-cut
            yAdd3 = obj.y(iin3).*0 - pi;
            yAdd = [yAdd12;yAdd3];
            iin = [iin12;iin3];
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.x > (pi/2+eps) | obj.x < (-pi/2-eps));
            signPos = sign(obj.x(outOfRangeInd));
            signPos(signPos == 0) = 1;
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - signPos.*pi;
            signPos = sign(obj.y(outOfRangeInd));
            signPos(signPos == 0) = 1;
            obj.y(outOfRangeInd) = signPos.*pi - obj.y(outOfRangeInd);
        end
    elseif type == 180
        if strcmp(obj.xRangeType,'pos')
            iout = find(abs(obj.y - 2*pi) < eps | (abs(obj.x - 0) < eps & obj.y > (pi/2+eps) & obj.y <= (3*pi/2+eps)));   % Redundant
            iin1 = find(abs(obj.y - 3*pi/2) < eps); 
            iin2 = find(abs(obj.y - pi/2) < eps & obj.x > eps); 
            iin = [iin1;iin2];
            xAdd = [obj.x(iin1);obj.x(iin2) + pi];
            yAdd = [obj.y(iin1) - 2*pi; pi - obj.y(iin2)];
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd1 = find(obj.y > (3*pi/2 + eps));
            obj.x(outOfRangeInd1) = obj.x(outOfRangeInd1);
            obj.y(outOfRangeInd1) = obj.y(outOfRangeInd1) - 2*pi;
            outOfRangeInd2 = find(obj.y > (pi/2 + eps));
            obj.x(outOfRangeInd2) = obj.x(outOfRangeInd2) + pi;
            obj.y(outOfRangeInd2) = pi - obj.y(outOfRangeInd2);
        elseif strcmp(obj.xRangeType,'sym')
%             iout = unique(find(abs(obj.y - pi) < eps | (abs(abs(obj.x) - pi/2) < eps & abs(obj.y) >= (pi/2-eps))));   % Redundant
            iout = find(abs(obj.y - pi) < eps | (abs(abs(obj.x) - pi/2) < eps & abs(obj.y) >= (pi/2-eps)));   % Redundant
            iin1 = find(abs(obj.x - 0) < eps & (obj.y >= (pi/2 - eps) | (obj.y <= -(pi/2 - eps) & obj.y > -(pi-eps))));
            iin2 = find(abs(abs(obj.y) - pi/2) < eps);
            iin = [iin1;iin2];
            signPos1 = sign(obj.x(iin1));
            signPos1(signPos1 == 0) = -1; % Add these to the x = 180 cut...
            xAdd1 = obj.x(iin1) - signPos1.*pi;
            signPos2 = sign(obj.x(iin2));
            signPos2(signPos2 == 0) = 1; 
            xAdd2 = obj.x(iin2) - signPos2.*pi;
            xAdd = [xAdd1;xAdd2];
            yAdd = sign(obj.y(iin)).*pi - obj.y(iin);
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.y > (pi/2 + eps) | obj.y < (-pi/2 - eps));
            signPos = sign(obj.x(outOfRangeInd));
            signPos(signPos == 0) = 1;
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - signPos.*pi;
            obj.y(outOfRangeInd) = sign(obj.y(outOfRangeInd)).*pi - obj.y(outOfRangeInd); % Never 0 so no need to worry about fixing that
        end
    end
else
    warning(['Cant shift a polar grid like ', obj.gridType, ' on a cartesian grid']);
end
% Sort
obj = obj.sortGrid;
% Update the object descriptive parameters
obj = setRangeTypes(obj);
end


function objNew = shiftRedun(obj,iout,iin,xAdd,yAdd)
objNew = obj;
Nredun = min(numel(iout),numel(iin));    % How many we have to remove/add

% Remove old redundant points
if Nredun > 0
    objNew.x(iout(1:Nredun)) = [];
    objNew.y(iout(1:Nredun)) = [];
    objNew.E1(iout(1:Nredun),:) = [];
    objNew.E2(iout(1:Nredun),:) = [];
    objNew.E3(iout(1:Nredun),:) = [];
end

% Add in new redundant points
E1Add = obj.E1(iin,:);
E2Add = obj.E2(iin,:);
E3Add = obj.E3(iin,:);
objNew.x = [objNew.x;xAdd(1:Nredun)];
objNew.y = [objNew.y;yAdd(1:Nredun)];
objNew.E1 = [objNew.E1;E1Add(1:Nredun,:)];
objNew.E2 = [objNew.E2;E2Add(1:Nredun,:)];
objNew.E3 = [objNew.E3;E3Add(1:Nredun,:)];
end