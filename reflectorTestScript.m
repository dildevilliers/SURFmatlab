% Script to scratch around while developing the code

clear all
close all

%% changeBase test
% Make global coor object
GC = coordinateSystem();
% GC.plot
% make a paraboloid
MR = reflector;
% Shift/rotate it
rotTh = deg2rad(-70);
rotPh = deg2rad(115);
rotPs = deg2rad(-45);
transVect = [0.5,0.3,-0.1];
MR.coor = MR.coor.rotGRASP([rotTh,rotPh,rotPs]);
MR.coor = MR.coor.translate(transVect);
MR.plot
MR.coor.plot
% Place feed at random position - in global coordinate system
feedCoor = coordinateSystem([1;0;0]);
feedCoor = feedCoor.rotGRASP([deg2rad(-135),0,0]);
feedCoor.plot
feedCoor1 = coordinateSystem([1;1;0]);
feedCoor1 = feedCoor1.rotGRASP([deg2rad(75),deg2rad(30),deg2rad(-76)]);
feedCoor1.plot


% % Get the required angles to rotate back to the GC
% Np = cross(MR.coor.z_axis,GC.z_axis);
% Np = Np./norm(Np);
% % plot3([0,Np(1)],[0,Np(2)],[0,Np(3)],'c')
% N = cross(Np,MR.coor.z_axis);
% N = N./norm(N);
% % plot3([0,N(1)],[0,N(2)],[0,N(3)],'m')
% % th = acos(dot(GC.z_axis,MR.coor.z_axis));
% th = angBetweenVectors(MR.coor.z_axis,GC.z_axis);
% ph = angBetweenVectors(MR.coor.x_axis,N);
% % Sort out sign
% phSign = sign(dot(cross(MR.coor.x_axis,N),MR.coor.z_axis));
% ph = ph*phSign;
% MRcoor = MR.coor;
% xyz_prime = MRcoor.rotGRASP(th,ph,0);
% ps = angBetweenVectors(xyz_prime.x_axis,GC.x_axis);
% % Sort out sign
% psSign = sign(dot(cross(xyz_prime.x_axis,GC.x_axis),GC.z_axis));
% ps = ps*psSign;
% % MRcoor = MRcoor.rotGRASP(th,ph,ps);
% % MRcoor.plot
[angGRASP] = getGRASPangBetweenCoors(GC,MR.coor);

% Shift/rotate back
figure
% change base of the start and end values of the coordinate vectors...
% fC2 = coordinateSystem(newPoints(:,1),newPoints(:,2)-newPoints(:,1),newPoints(:,3)-newPoints(:,1));
fC2 = changeBase([feedCoor,feedCoor1],GC,MR.coor);
fC2(1).plot
fC2(2).plot

MR.coor = MR.coor.translate(-transVect);
MR.coor = rotGRASP(MR.coor,angGRASP);
MR.plot
MR.coor.plot

