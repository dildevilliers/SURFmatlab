%Name: testScript_FarField.m
%Description:
%   Script to test the Farfield object. A default Farfield is created and the
%   use cases of its various methods are shown.

% A default FarField class :
% Generates a Gaussian beam pattern at a single frequency (1 GHz) over the full sphere, 5 degree angular resolution

close all

% FF = FarField;
FF = FarField.readFEKOffe([pwd, '\testScripts\patch']);

%FF = FF.pol2circular();
% figure 
% FF.plot;

gridType = 'AzEl';
coorType = 'Ludwig3';
polType = 'circular';
xRangeType = 'pos';

handleGridType = str2func(['grid2',gridType]);
handleCoorType = str2func(['coor2',coorType]);
handlePolType = str2func(['pol2',polType]);
FF = handleGridType(FF);
FF = handleCoorType(FF,0);
FF = handlePolType(FF);
FF = FF.setXrange(xRangeType);
figure
%FF.plot('plotType','2D','step',1,'showGrid',true,'output','E1','outputType','mag');
% figure
% FF.plot('plotType','2D','step',1,'showGrid',true,'output','E2','outputType','mag');

% figure
% plotType = 'polar'; % '3D', '2D', 'polar', 'cartesian'
% FF.plot('plotType', plotType);

% figure
% FF.plotPrincipleCuts;
% 
 figure 
 %FF.plot('plotType', 'polar', 'output', 'CO_XP', 'outputType', 'mag')

FF2 = grid2PhTh(FF);
figure
plotGrid(FF2);
FF.plot();
% 
% FF3 = grid2AzEl(FF);
% figure
% plotGrid(FF3);
% 
% FF4 = grid2ArcSin(FF);
% figure
% plotGrid(FF4);
% 
% FF5 = grid2ElAz(FF);
% figure
% plotGrid(FF5);
% 
% FF6 = grid2TrueView(FF);
% figure
% plotGrid(FF6);


