function writeGaussianBeamInfo2TOR(fileID,GaussName,frequencyVal,coor_sysVal,taper_angleVal,taperVal)
fprintf(fileID, '%s  gaussian_beam_pattern  \n',GaussName);
fprintf(fileID, '(\n');
writeTORstandardLine(fileID,'frequency',frequencyVal,0,1,0);
writeTORstandardLine(fileID,'coor_sys',coor_sysVal,0,1,0);
writeTORstandardLine(fileID,'taper_angle',taper_angleVal,0,1,1);
writeTORstandardLine(fileID,'taper',taperVal,0,0,1);
fprintf(fileID, ')\n\n');
