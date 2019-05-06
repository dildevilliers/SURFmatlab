function exportSingleReflectorToTOR(obj,fullpathName,freqValGhz,prefixName)
fileID = fopen(fullpathName,'wt');

frequencyVal = freqValGhz;
prefixString = prefixName;

single_reflectorCoorAngGRASP = obj.PR.coor.getGRASPangles();
single_reflectorFEEDCoorAngGRASP = obj.feedCoor.getGRASPangles();

reflectorSurfaceName = strcat(prefixString,'surface');
reflectorRimName = strcat(prefixString,'rim');
reflectorCoorName = strcat(prefixString,'global_coor');
reflectorName = strcat(prefixString,'reflector');
feedCoorName = strcat(prefixString,'feed_coor');
frequencyName = strcat(prefixString,'frequencies');
feedName = strcat(prefixString,'feed');
POName = strcat(prefixString,'po');
CutName = strcat(prefixString,'cut');
CutCoorName = strcat(prefixString,'cut_coor');

writeCoorSysGraspAngles2TOR(fileID,reflectorCoorName, obj.PR.coor.origin.x,obj.PR.coor.origin.y,obj.PR.coor.origin.z,rad2deg(single_reflectorCoorAngGRASP(1)),rad2deg(single_reflectorCoorAngGRASP(2)),rad2deg(single_reflectorCoorAngGRASP(3)),obj.PR.coor.base);
%=========================================================================
%HARD CODED FREQ. NO FREQUENCY INFO GIVEN IN GEOMETRY CLASS. BELOW.
%=========================================================================
writeFreqInfo2TOR(fileID,frequencyName,frequencyVal);
%=========================================================================
%HARD CODED FREQ. NO FREQUENCY INFO GIVEN IN GEOMETRY CLASS. ABOVE.
%=========================================================================
writeSurfaceParaboloid2TOR(fileID,reflectorSurfaceName, obj.PR.surface.focalLength,obj.PR.surface.vertex.x,obj.PR.surface.vertex.y,obj.PR.surface.vertex.z);
writeRimElliptical2TOR(fileID,reflectorRimName,obj.PR.rim.centre(1),obj.PR.rim.centre(2),obj.PR.rim.halfAxis(1),obj.PR.rim.halfAxis(2),obj.PR.rim.rotation);
writeReflector2TOR(fileID,reflectorName, reflectorCoorName, reflectorSurfaceName, reflectorRimName);
writeCoorSysGraspAngles2TOR(fileID,feedCoorName, obj.feedCoor.origin.x,obj.feedCoor.origin.y,obj.feedCoor.origin.z,rad2deg(single_reflectorFEEDCoorAngGRASP(1)),rad2deg(single_reflectorFEEDCoorAngGRASP(2)),rad2deg(single_reflectorFEEDCoorAngGRASP(3)),reflectorCoorName);
%=========================================================================
%START OF ELECTRICAL INFORMATION SECTION. NO INFO GIVEN IN GEOMETRY CLASS
%REGARDING THIS, SO THIS IS HARDCODED. BELOW.
%=========================================================================
writeGaussianBeamInfo2TOR(fileID,feedName,frequencyName,feedCoorName,45.2397302805424,-12.0);
writeCoorSysGraspAngles2TOR(fileID,CutCoorName, 0,0,0,0,0,0,reflectorCoorName);
writeCutInfo2TOR(fileID,CutName,CutCoorName,-7.15701779967802,7.15701779967802,161,0,90,3,frequencyName);
writePOInfo2TOR(fileID,POName,frequencyName,reflectorName,reflectorCoorName);
%=========================================================================
%END OF ELECTRICAL INFORMATION SECTION. NO INFO GIVEN IN GEOMETRY CLASS
%REGARDING THIS, SO THIS IS HARDCODED. ABOVE.
%=========================================================================
%=========================================================================
%THIS NEXT SECTION HARDCODES THE GRAPHICS INFO. THIS IS BECAUSE THERE'S NO
%DYNAMIC ADJUSTMENT OF IT IN THE GEOMETRY CLASS. NO INFORMATION ABOUT IT
%CAN BE CREATED, SO IT CAN BE HARDCODED. HAPPENS BELOW.
%=========================================================================
writeStaticBasicGraphicsViewInfo2TOR(fileID);
%=========================================================================
%THIS NEXT SECTION HARDCODES THE GRAPHICS INFO. THIS IS BECAUSE THERE'S NO
%DYNAMIC ADJUSTMENT OF IT IN THE GEOMETRY CLASS. NO INFORMATION ABOUT IT
%CAN BE CREATED, SO IT CAN BE HARDCODED. HAPPENS ABOVE.
%=========================================================================
fclose(fileID);