function writePOInfo2TOR(fileID,POName,frequencyVal,scattererVal,coor_sysVal)
fprintf(fileID, '%s  po_single_face_scatterer  \n',POName);
fprintf(fileID, '(\n');
writeTORstandardLine(fileID,'frequency',frequencyVal,0,1,0);
writeTORstandardLine(fileID,'scatterer',scattererVal,0,1,0);
writeTORstandardLine(fileID,'coor_sys',coor_sysVal,0,0,0);
fprintf(fileID, ')\n\n');