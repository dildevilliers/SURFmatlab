function writeRimElliptical2TOR(fileID,RimName,centreX,centreY,radiusHA_x,radiusHA_y,rotation)
fprintf(fileID, '%s  elliptical_rim  \n',RimName);
fprintf(fileID, '(\n');
writeTOR2PartLine(fileID,'centre','x',centreX,'y',centreY,1,1);
writeTOR2PartLine(fileID,'half_axis','x',radiusHA_x,'y',radiusHA_y,1,1);
writeTORstandardLine(fileID,'rotation',rotation,0,0);
fprintf(fileID, ')\n\n');