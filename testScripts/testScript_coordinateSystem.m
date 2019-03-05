% Script to test the coordinateSystem class

close all
clear all

originBase = [-2,0.6,1].';
angGRASPBase = deg2rad([45,30,0]);
% originBase = [0,0,0].';
% angGRASPBase = deg2rad([0,0,0].');

originCoor = [-1,3,2];
angGRASPCoor = deg2rad([-30,67,80]);
% originCoor = [2,2,0];
% angGRASPCoor = deg2rad([0,0,90]);

originOtherBase = [2,-1,-0.5].';
angGRASPotherBase = deg2rad([20,-78,143]);

coorGlob = coordinateSystem();
coorGlob.plot
hold on

coorBase = coordinateSystem(originBase);
coorBase = coorBase.rotGRASP(angGRASPBase);
coorBase.plot

coorSys = coordinateSystem(originCoor);
coorSys = coorSys.rotGRASP(angGRASPCoor);
coorSys.base = coorBase;
coorSys.plot

coorSysInGlobal = coorSys.getInGlobal;
coorSysInGlobal.origin;
coorSysInGlobal.x_axis;
coorSysInGlobal.y_axis;

coorOtherBase = coordinateSystem(originOtherBase);
coorOtherBase = coorOtherBase.rotGRASP(angGRASPotherBase);
coorSysInOtherBase = coorSys.redefineToOtherBase(coorOtherBase);
coorSysInOtherBase.origin
coorSysInOtherBase.x_axis
coorSysInOtherBase.y_axis


% U = [1,1,0].';
% Uprime = changeBase(U,coorSys,coorBase)
