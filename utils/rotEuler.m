function Xd = rotEuler(X,angEuler)

% function Xd = rotEuler(X,angEuler)
% Returns the Euler angle (angEuler) rotated vector Xd = [x';y';z'] of the input
% vector X = [x;y;z] 
% Vectors can be rows of equal length

coor_base = coordinateSystem();     % Work from global coordinate system
% Make a GRASP rotated coordinate system
coor_new = coordinateSystem();
coor_new = rotEuler(coor_new,angEuler);
Xd = changeBase(X,coor_new,coor_base);
