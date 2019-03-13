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
[FFM,MaskPointing,M] = SP.getMask;
figure
plot(FFM,'output','Directivity','outputType','mag','plotType','2D','scaleMag','lin','step',1,'norm',1)
% Also plot the mask pointing
[x,y,z] = deal(zeros(size(MaskPointing)));
for mm = 1:1:length(MaskPointing)
    endPoints = MaskPointing(mm).origin.addVect(MaskPointing(mm).z_axis);
    x(mm) = endPoints.x;
    y(mm) = endPoints.y;
    z(mm) = endPoints.z;
end
figure
plot3(x(M),y(M),z(M),'k.'), grid on, hold on
plot3(x(~M),y(~M),z(~M),'r.')
axis equal