close all
clear all

SP = symmetricParaboloid;
% SP.plot
% figure
% SP.plot3D
% figure
% SP.plotRayTrace
% pathLengthStruct = SP.getPathLength
FFM = SP.getMask;
plot(FFM,'output','Directivity','outputType','mag','plotType','2D','scaleMag','lin','step',1,'norm',1)