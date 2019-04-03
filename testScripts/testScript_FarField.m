%Name: testScript_FarField.m
%Description:
%   Script to test the Farfield object. A patch antenna Farfield is read in from a FEKO ffe file
%   and the use cases of its various methods are shown.

close all
FF = FarField.readFEKOffe([pwd, '\testScripts\patch']);

gridType = {'PhTh' 'DirCos' 'AzEl' 'ElAz' 'TrueView' 'ArcSin'}; 
coorType = {'spherical' 'Ludwig1' 'Ludwig2AE' 'Ludwig2EA' 'Ludwig3'};
polType = {'linear' 'circular' 'slant'};
xRangeType = 'pos';

for i = 1:6
    handleGridType = str2func(['grid2',gridType{i}]);
    FF = handleGridType(FF);
    figure
    plotGrid(FF);
end

for i = 1:6
    k = i;
    if i > 5
        k = 1;
    end
    for j = 1:3
        handleGridType = str2func(['grid2',gridType{i}]);
        handleCoorType = str2func(['coor2',coorType{k}]);
        handlePolType = str2func(['pol2',polType{j}]);
        FF = handleGridType(FF);
        FF = handleCoorType(FF,0);
        FF = handlePolType(FF);
        FF = FF.setXrange(xRangeType);
        figure
        FF.plot('plotType','2D','step',1,'showGrid',true,'output','E1','outputType','mag');
        figure
        FF.plot('plotType','2D','step',1,'showGrid',true,'output','E2','outputType','mag');
%         figure
%         FF.plot('plotType','2D','step',1,'showGrid',true,'output','CO_XP','outputType','mag');
%         figure
%         FF.plot('plotType','2D','step',1,'showGrid',true,'output','XP_CO','outputType','mag');
    end  
end

figure
plotType = 'polar'; % '3D', '2D', 'polar', 'cartesian'
FF.plot('plotType', plotType);

figure
FF.plotPrincipleCuts;

