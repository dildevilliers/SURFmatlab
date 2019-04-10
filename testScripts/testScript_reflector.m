% Script to scratch around while developing the code
clear all
close all


%% Basic ellipsoid test
GC = coordinateSystem;
SR = reflector;
SR.surface = ellipsoid(5,4);
SR.plot
SR.coor.plot

keyboard
%% changeBase test

figure
% Make global coor object
GC = coordinateSystem();
% GC.plot
% make a reflector
MR = reflector;
% % Test ellipsoid
% MR.surface = ellipsoid(5,4);
% Test hyperboloid
MR.surface = hyperboloid(4,5);
% Shift/rotate it
rotTh = deg2rad(45);
rotPh = deg2rad(-123);
rotPs = deg2rad(65);
transVect = [0.5,0.3,-0.1].';
MR.coor = MR.coor.rotGRASP([rotTh,rotPh,rotPs]);
MR.coor = MR.coor.translate(transVect);
MR.plot(10000,'cart')
MR.coor.plot
% Place feed at random position - in global coordinate system
feedCoor = coordinateSystem(pnt3D(1,0,0));
feedCoor = feedCoor.rotGRASP([deg2rad(-135),0,0]);
feedCoor.plot
feedCoor1 = coordinateSystem(pnt3D(1,1,1));
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
close all
clear all
% S = paraboloid(pnt3D(0,0,0),5);
% R = ellipticalRim([1;2],[0.5;1]);
% S = ellipsoid(5,4);
% R = ellipticalRim([0;0],[1.5;1.5]);
S = hyperboloid(4,5);
R = ellipticalRim([0;0],[1.5;1.5]);
C = coordinateSystem();
C = C.rotGRASP(deg2rad([0,0,0]));
R = reflector(S,R,C);
R.plot(10000,'polarThin')
% R.plot(1000,'x0')
% R.plot(1000,'x0')
R.plotNorms(1000,100,1,'polarThin')

%% Masking/ray tracing test

close all
clear all

S = paraboloid(pnt3D(0,0,0),5);
R = ellipticalRim([0;0],[5;5]);
C = coordinateSystem(pnt3D(0,0,5));
C = C.rotGRASP(deg2rad([-90,0,0]));
C0 = coordinateSystem();
C0 = C0.rotGRASP(deg2rad([0,0,0]));

% S = paraboloid(pnt3D(2,2,0),5);
% R = ellipticalRim([2;2],[5;4]);
% C = coordinateSystem(pnt3D(1,3,5));
% C = C.rotGRASP(deg2rad([-145,30,90]));
% C0 = coordinateSystem();
% C0 = C0.rotGRASP(deg2rad([60,-45,0]));

% S = ellipsoid(5,4);
% R = ellipticalRim([0;0],[1.2;1.2]);
% C = coordinateSystem(pnt3D(0,0,0));
% C = C.rotGRASP(deg2rad([-15,45,0]));
% C0 = coordinateSystem();
% C0 = C0.rotGRASP(deg2rad([60,-45,0]));

% S = hyperboloid(3.5,5);
% R = ellipticalRim([0;0],[1.2;1.2]);
% C = coordinateSystem(pnt3D(0,0,-2*S.f));
% C = C.rotGRASP(deg2rad([-15,45,0]));
% C0 = coordinateSystem();
% C0 = C0.rotGRASP(deg2rad([60,-45,0]));

R = reflector(S,R,C0);
R.plot
C.plot
% C0.plot
% R.getMaskFunction(C);
th = linspace(-pi,pi,200);
ph = ones(size(th)).*deg2rad(45);
FF = FarField;
M = R.getMask(C,[ph.',th.'],1);
% R.plotMask(C,[ph.',th.'],12)
% R.plotMask(C,FF,12)

% R.getRayInterceptPoint(C,ph.',th.',500)
% Pintercept = R.getRayInterceptPoint(C,FF.ph,FF.th,500);
R.reflectRays(C,ph.',th.',200,1)
