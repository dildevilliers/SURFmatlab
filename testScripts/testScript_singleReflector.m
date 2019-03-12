close all
clear all

% SP = symmetricParaboloid;
SP = singleReflector(1,1.7,2.5);
SP.plot
% figure
% SP.plot3D
figure
SP.plotRayTrace
% pathLengthStruct = SP.getPathLength
FFM = SP.getMask;
figure
plot(FFM,'output','Directivity','outputType','mag','plotType','2D','scaleMag','lin','step',1,'norm',1)