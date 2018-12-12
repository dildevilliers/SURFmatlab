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

if strcmp(obj.gridType,'PhTh')
    if type == 360
        if strcmp(obj.xRangeType,'pos')
            % Redundant pole: remove ph > 180, th = 180; add ph = 0, th = 180:360
            iout = find(obj.x  > (pi+eps) & abs(obj.y - pi) < 1e-10);   % Redundant
            iin = find(abs(obj.x - 0) < 1e-10);                   % Will become redundant after inserting
            xAdd = obj.x(iin);
            yAdd = 2*pi - obj.y(iin);
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.x > (pi+eps));
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - pi;
            obj.y(outOfRangeInd) = 2*pi - obj.y(outOfRangeInd);
        elseif strcmp(obj.xRangeType,'sym')
            % Redundant pole: remove |ph| > 90, th = 0 and  ph = -180, th = 0:180; add ph = abs(90), th = 0:180
            iout = unique([find((obj.x > (pi/2+eps) | obj.x < (-pi/2-eps)) & abs(obj.y - 0) < 1e-10); find(abs(obj.x + pi) < 1e-10)]);   % Redundant
            iin1 = find(abs(obj.x + pi/2) < 1e-10);
            iin2 = find(abs(obj.x - pi/2) < 1e-10);
            iin = [iin1;iin2];
            xAdd1 = obj.x(iin1);
            xAdd2 = obj.x(iin2);
            xAdd = [xAdd1;xAdd2];
            yAdd = -obj.y(iin);
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.x > (pi/2+eps) | obj.x < (-pi/2-eps));
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - sign(obj.x(outOfRangeInd)).*pi;
            obj.y(outOfRangeInd) = -obj.y(outOfRangeInd);
        end
    elseif type == 180
        if strcmp(obj.xRangeType,'pos')
            % Redundant pole: remove ph = 180, th = 0:180; add ph = 0:180, th = 180
            iout = find(obj.y < (pi-eps) & abs(obj.x - pi) < 1e-10);   % Redundant
            iin = find(obj.x < (pi-eps) & abs(obj.y - pi) < 1e-10);    % Will become redundant after inserting
            xAdd = obj.x(iin);
            yAdd = obj.y(iin);
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.y > (pi+eps));
            obj.y(outOfRangeInd) = 2*pi - obj.y(outOfRangeInd);
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) + pi;
        elseif strcmp(obj.xRangeType,'sym')
            % Redundant pole: remove ph = |90|, th < 0; add ph > |90|, th = 0
            iout = unique(find((abs(obj.x - pi/2) < 1e-10 | abs(obj.x + pi/2) < 1e-10) & (obj.y < (pi/2-eps))));   % Redundant
            iin1 = find(obj.x < -eps & abs(obj.y - 0) < 1e-10);
            iin2 = find(obj.x > eps & abs(obj.y - 0) < 1e-10);
            iin = [iin1;iin2];
            xAdd1 = -pi - obj.x(iin1);
            xAdd2 = pi - obj.x(iin2);
            xAdd = [xAdd1;xAdd2];
            yAdd = obj.y(iin);
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.y < -eps);
            obj.y(outOfRangeInd) = -obj.y(outOfRangeInd);
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - sign(obj.x(outOfRangeInd)).*pi;
        end
    end
elseif strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
    % ToDo
    if type == 360
        if strcmp(obj.xRangeType,'pos')
            % Redundant pole: remove ph > 180, th = 180; add ph = 0, th = 180:360
            iout = find(obj.x  > (pi+eps) & abs(obj.y - pi) < 1e-10);   % Redundant
            iin = find(abs(obj.x - 0) < 1e-10);                   % Will become redundant after inserting
            xAdd = obj.x(iin);
            yAdd = 2*pi - obj.y(iin);
            obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
            % Apply shift
            outOfRangeInd = find(obj.x > (pi+eps));
            obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - pi;
            obj.y(outOfRangeInd) = 2*pi - obj.y(outOfRangeInd);
        elseif strcmp(obj.xRangeType,'sym')
%             % Redundant pole: remove |ph| > 90, th = 0 and  ph = -180, th = 0:180; add ph = abs(90), th = 0:180
%             iout = unique([find((obj.x > (pi/2+eps) | obj.x < (-pi/2-eps)) & abs(obj.y - 0) < 1e-10); find(abs(obj.x + pi) < 1e-10)]);   % Redundant
%             iin1 = find(abs(obj.x + pi/2) < 1e-10);
%             iin2 = find(abs(obj.x - pi/2) < 1e-10);
%             iin = [iin1;iin2];
%             xAdd1 = obj.x(iin1);
%             xAdd2 = obj.x(iin2);
%             xAdd = [xAdd1;xAdd2];
%             yAdd = -obj.y(iin);
%             obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
%             % Apply shift
%             outOfRangeInd = find(obj.x > (pi/2+eps) | obj.x < (-pi/2-eps));
%             obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - sign(obj.x(outOfRangeInd)).*pi;
%             obj.y(outOfRangeInd) = -obj.y(outOfRangeInd);
        end
    elseif type == 180
        if strcmp(obj.xRangeType,'pos')
%             % Redundant pole: remove ph = 180, th = 0:180; add ph = 0:180, th = 180
%             iout = find(obj.y < (pi-eps) & abs(obj.x - pi) < 1e-10);   % Redundant
%             iin = find(obj.x < (pi-eps) & abs(obj.y - pi) < 1e-10);    % Will become redundant after inserting
%             xAdd = obj.x(iin);
%             yAdd = obj.y(iin);
%             obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
%             % Apply shift
%             outOfRangeInd = find(obj.y > (pi+eps));
%             obj.y(outOfRangeInd) = 2*pi - obj.y(outOfRangeInd);
%             obj.x(outOfRangeInd) = obj.x(outOfRangeInd) + pi;
        elseif strcmp(obj.xRangeType,'sym')
%             % Redundant pole: remove ph = |90|, th < 0; add ph > |90|, th = 0
%             iout = unique(find((abs(obj.x - pi/2) < 1e-10 | abs(obj.x + pi/2) < 1e-10) & (obj.y < (pi/2-eps))));   % Redundant
%             iin1 = find(obj.x < -eps & abs(obj.y - 0) < 1e-10);
%             iin2 = find(obj.x > eps & abs(obj.y - 0) < 1e-10);
%             iin = [iin1;iin2];
%             xAdd1 = -pi - obj.x(iin1);
%             xAdd2 = pi - obj.x(iin2);
%             xAdd = [xAdd1;xAdd2];
%             yAdd = obj.y(iin);
%             obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
%             % Apply shift
%             outOfRangeInd = find(obj.y < -eps);
%             obj.y(outOfRangeInd) = -obj.y(outOfRangeInd);
%             obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - sign(obj.x(outOfRangeInd)).*pi;
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


% % Function to remove the redundant values
% function [objSmall,Nredun] = removeRedun(objFull,iout,iin)
% objSmall = objFull;
% Nredun = min(numel(iout),numel(iin));    % How many we have to remove/add
% % Truncate the indexes to the shortest length
% iout = iout(1:Nredun);
% % iinSmall = iin(1:Nredun);
% if Nredun > 0
%     objSmall.x(iout) = [];
%     objSmall.y(iout) = [];
%     objSmall.E1(iout,:) = [];
%     objSmall.E2(iout,:) = [];
%     objSmall.E3(iout,:) = [];
% end
% end

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