function exportDualReflectorToTOR(obj,fullpathName,freqValGhz,prefixName)
fileID = fopen(fullpathName,'wt');
%{
figure
obj.plot3D(10000,[1,1,1])
figure
obj.plot
figure
obj.plotRayTrace
%}
frequencyVal = freqValGhz;
prefixString = prefixName;

dual_reflectorPRCoorAngGRASP = obj.PR.coor.getGRASPangles();
dual_reflectorSRCoorAngGRASP = obj.SR.coor.getGRASPangles();
dual_reflectorFEEDCoorAngGRASP = obj.feedCoor.getGRASPangles();

reflectorCoorName = strcat(prefixString,'global_coor');
reflectorSurfaceName = strcat(prefixString,'main_surface');
reflectorRimName = strcat(prefixString,'main_rim');
reflectorName = strcat(prefixString,'main_reflector');
frequencyName = strcat(prefixString,'frequencies');
feedCoorName = strcat(prefixString,'feed_coor');
SUBreflectorCoorName = strcat(prefixString,'sub_coor');
SUBreflectorSurfaceName = strcat(prefixString,'sub_surface');
SUBreflectorRimName = strcat(prefixString,'sub_rim');
SUBreflectorName = strcat(prefixString,'sub_reflector');

feedName = strcat(prefixString,'feed');
POName = strcat(prefixString,'main_po');
SUBPOName = strcat(prefixString,'sub_po');
CutName = strcat(prefixString,'cut');
CutCoorName = strcat(prefixString,'cut_coor');

writeCoorSysGraspAngles2TOR(fileID,reflectorCoorName, obj.PR.coor.origin.x,obj.PR.coor.origin.y,obj.PR.coor.origin.z,rad2deg(dual_reflectorPRCoorAngGRASP(1)),rad2deg(dual_reflectorPRCoorAngGRASP(2)),rad2deg(dual_reflectorPRCoorAngGRASP(3)),obj.PR.coor.base);
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
writeCoorSysGraspAngles2TOR(fileID,feedCoorName, obj.feedCoor.origin.x,obj.feedCoor.origin.y,obj.feedCoor.origin.z,rad2deg(dual_reflectorFEEDCoorAngGRASP(1)),rad2deg(dual_reflectorFEEDCoorAngGRASP(2)),rad2deg(dual_reflectorFEEDCoorAngGRASP(3)),reflectorCoorName);
if obj.sigma == -1 %CASSEGRAIN SYSTEM.
    writeSurfaceHyperboloid2TOR(fileID,SUBreflectorSurfaceName, obj.SR.surface.vertexDistance,obj.SR.surface.fociDistance);
elseif obj.sigma == 1 %GREGORIAN SYSTEM.
    writeSurfaceEllipsoid2TOR(fileID,SUBreflectorSurfaceName, obj.SR.surface.vertexDistance,obj.SR.surface.fociDistance,obj.SR.surface.rotAng);
end
writeRimElliptical2TOR(fileID,SUBreflectorRimName,obj.SR.rim.centre(1),obj.SR.rim.centre(2),obj.SR.rim.halfAxis(1),obj.SR.rim.halfAxis(2),obj.SR.rim.rotation);
writeReflector2TOR(fileID,SUBreflectorName, SUBreflectorCoorName, SUBreflectorSurfaceName, SUBreflectorRimName);
if obj.sigma == -1 %CASSEGRAIN SYSTEM.
    writeCoorSysGraspAngles2TOR(fileID,SUBreflectorCoorName, obj.SR.coor.origin.x,obj.SR.coor.origin.y,obj.SR.coor.origin.z,rad2deg(dual_reflectorSRCoorAngGRASP(1)+pi),rad2deg(dual_reflectorSRCoorAngGRASP(2)),rad2deg(dual_reflectorSRCoorAngGRASP(3)+pi),obj.SR.coor.base);
elseif obj.sigma == 1 %GREGORIAN SYSTEM.
    writeCoorSysGraspAngles2TOR(fileID,SUBreflectorCoorName, obj.SR.coor.origin.x,obj.SR.coor.origin.y,obj.SR.coor.origin.z,rad2deg(dual_reflectorSRCoorAngGRASP(1)),rad2deg(dual_reflectorSRCoorAngGRASP(2)),rad2deg(dual_reflectorSRCoorAngGRASP(3)),obj.SR.coor.base);
end
%=========================================================================
%START OF ELECTRICAL INFORMATION SECTION. NO INFO GIVEN IN GEOMETRY CLASS
%REGARDING THIS, SO THIS IS HARDCODED. BELOW.
%=========================================================================
writeGaussianBeamInfo2TOR(fileID,feedName,frequencyName,feedCoorName,10.6884202229611,-12.0);
writeCoorSysGraspAngles2TOR(fileID,CutCoorName, 0,0,0,0,0,0,reflectorCoorName);
writeCutInfo2TOR(fileID,CutName,CutCoorName,-2.36921968541065,2.36921968541065,161,0,90,3,frequencyName);
writePOInfo2TOR(fileID,POName,frequencyName,reflectorName,reflectorCoorName);
writePOInfo2TOR(fileID,SUBPOName,frequencyName,SUBreflectorName,reflectorCoorName);
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