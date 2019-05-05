function writeSurfaceParaboloid2TOR(fileID,SurfaceParaboloidName, SurfaceParaboloidFocalLength,vertexX,vertexY,vertexZ)
fprintf(fileID, '%s  paraboloid  \n',SurfaceParaboloidName);
fprintf(fileID, '(\n');
writeTORstandardLine(fileID,'focal_length',SurfaceParaboloidFocalLength,1,1);
writeTOR3PartLine(fileID,'vertex','x',vertexX,'y',vertexY,'z',vertexZ,1,0);
fprintf(fileID, ')\n\n');