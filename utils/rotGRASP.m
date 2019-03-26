function Xd = rotGRASP(X,angGRASP)

% function Xd = rotGRASP(X,angGRASP)
% Returns the GRASP angle (angGRASP) rotated vector Xd = [x';y';z'] of the input
% vector X = [x;y;z] 
% Vectors can be rows of equal length

coor_base = coordinateSystem();     % Work from global coordinate system
% Make a GRASP rotated coordinate system
coor_new = coordinateSystem();
coor_new = rotGRASP(coor_new,angGRASP);
Xd = changeBase(X,coor_new,coor_base);
