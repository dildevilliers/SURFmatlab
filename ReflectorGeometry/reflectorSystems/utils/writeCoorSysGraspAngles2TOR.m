function writeCoorSysGraspAngles2TOR(fileID,CoorSysName, originX,originY,originZ,theta,phi,psi,base)
fprintf(fileID, '%s  coor_sys_grasp_angles  \n',CoorSysName);
fprintf(fileID, '(\n');
writeTOR3PartLine(fileID,'origin','x',originX,'y',originY,'z',originZ,1,1);
writeTORstandardLine(fileID,'theta',theta,0,1);
writeTORstandardLine(fileID,'phi',phi,0,1);
if ~isempty(base)
    writeTORstandardLine(fileID,'psi',psi,0,1);
    writeTORstandardLine(fileID,'base',base,0,0,0);
else
    writeTORstandardLine(fileID,'psi',psi,0,0);
end
fprintf(fileID, ')\n\n');