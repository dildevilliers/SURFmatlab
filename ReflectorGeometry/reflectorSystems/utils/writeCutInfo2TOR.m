function writeCutInfo2TOR(fileID,CutName,coor_sysVal,startVal1,endVal1,npVal1,startVal2,endVal2,npVal2,frequencyVal)
fprintf(fileID, '%s  spherical_cut  \n',CutName);
fprintf(fileID, '(\n');
writeTORstandardLine(fileID,'coor_sys',coor_sysVal,0,1,0);
writeTOR3PartLine(fileID,'theta_range','start',startVal1,'end',endVal1,'np',npVal1,0,1);
writeTOR3PartLine(fileID,'phi_range','start',startVal2,'end',endVal2,'np',npVal2,0,1);
fprintf(fileID, '  comment          : "Field data in cuts",\n');
writeTORstandardLine(fileID,'frequency',frequencyVal,0,0,0);
fprintf(fileID, ')\n\n');