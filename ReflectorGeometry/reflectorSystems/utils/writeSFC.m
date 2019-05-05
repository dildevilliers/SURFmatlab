function writeSFC(pathNameExt,informationName, numOfPoints,X,Y,Z)
% HOW IT WORKS [file_format](sfc for grasp):
%IRREGULAR XY-GRID, PSEUDO SPLINES:
%1) TEXT -> RECORD WITH IDENIFICATION TEXT.
%2) N_POINTS -> NUMBER OF (X,Y,Z) VALUES.
%3) X(I), Y(I), Z(I), I=1, N_POINTS (3 REAL NUMBERS).
%
%EXAMPLE:
%INTERP_Shape_Cut_Paraboloid_Data_
%2556
%0.059713 -0.105657 0.053975
%0.059713 -0.097389 0.052713
%...
fileID = fopen(pathNameExt,'wt');
%%THIS SECTION ATTEMPTS TO WRITE A .SFC FILE TO GRASP.
%=========================================================================
fprintf(fileID, '%s\n',informationName);
fprintf(fileID, '%d\n',numOfPoints);
for q=1:1:numOfPoints
    fprintf(fileID, '%f %f %f\n',X(q),Y(q),Z(q));
end
%=========================================================================
fclose(fileID);
