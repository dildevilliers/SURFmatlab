function writeSurfaceEllipsoid2TOR(fileID,SurfaceEllipsoidName, vertex_distanceVal,foci_distanceVal,axis_angleVal)
fprintf(fileID, '%s  ellipsoid  \n',SurfaceEllipsoidName);
fprintf(fileID, '(\n');
writeTORstandardLine(fileID,'vertex_distance',vertex_distanceVal,1,1);
writeTORstandardLine(fileID,'foci_distance',foci_distanceVal,1,1);
writeTORstandardLine(fileID,'axis_angle',axis_angleVal,0,0);
fprintf(fileID, ')\n\n');