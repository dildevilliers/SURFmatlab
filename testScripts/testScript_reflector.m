% Script to scratch around while developing the code

clear all
close all

%% changeBase test
figure
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
MR.plot(10000,'cart')
MR.coor.plot
% Place feed at random position - in global coordinate system
feedCoor = coordinateSystem([1;0;0]);
feedCoor = feedCoor.rotGRASP([deg2rad(-135),0,0]);
feedCoor.plot
feedCoor1 = coordinateSystem([1;1;1]);
feedCoor1 = feedCoor1.rotGRASP([deg2rad(75),deg2rad(30),deg2rad(-76)]);
feedCoor1.plot

% Shift/rotate back to global - by plotting the coordinate systems in the
% MR coordinate system
figure
% fC2 = changeBase([feedCoor,feedCoor1],MR.coor,GC);
fC2(1) = feedCoor.redefineToOtherBase(MR.coor);
fC2(2) = feedCoor1.redefineToOtherBase(MR.coor);
fC2(1).plotLocal
fC2(2).plotLocal
GC.plot

% % Get the required angles to rotate back to the GC
[angGRASP] = getGRASPangBetweenCoors(GC,MR.coor);
MR.coor = MR.coor.translate(-transVect);
MR.coor = rotGRASP(MR.coor,angGRASP);
MR.plot
MR.coor.plot

%% Grid and plotting test
R = reflector();
% R.plot(10000,'polarThin')
% R.plot(1000,'x0')
% R.plot(1000,'x0')
R.plotNorms(1000,200,'polarThin')

%% Masking test
S = paraboloid([0;0;-5],5);
R = ellipticalRim([0;0],[5;5]);
C0 = coordinateSystem();
C0 = C0.rotGRASP(deg2rad([0,0,0]));
R = reflector(S,R,C0);
C = coordinateSystem([2,6,0]);
C = C.rotGRASP(deg2rad([-145,25,90]));
R.plot
C.plot
R.getMaskFunction(C)
