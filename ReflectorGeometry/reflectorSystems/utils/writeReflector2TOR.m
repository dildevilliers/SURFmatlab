function writeReflector2TOR(fileID,reflectorName, coorSys, surface, rim)
fprintf(fileID, '%s  reflector  \n',reflectorName);
fprintf(fileID, '(\n');
writeTORstandardLine(fileID,'coor_sys',coorSys,0,1,0);

surface = parseThisString(surface);
fprintf(fileID, '  surfaces         : sequence(%s),\n', surface);

writeTORstandardLine(fileID,'rim',rim,0,0,0);
fprintf(fileID, ')\n\n');