function writeSurfaceHyperboloid2TOR(fileID,SurfaceHyperboloidName, vertex_distanceVal,foci_distanceVal)
fprintf(fileID, '%s  hyperboloid  \n',SurfaceHyperboloidName);
fprintf(fileID, '(\n');
writeTORstandardLine(fileID,'vertex_distance',vertex_distanceVal,1,1);
writeTORstandardLine(fileID,'foci_distance',foci_distanceVal,1,0);
fprintf(fileID, ')\n\n');