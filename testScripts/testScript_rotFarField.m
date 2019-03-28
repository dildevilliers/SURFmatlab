% Script to test FarField rotations
clear all
close all

plotDim = 3;    % Select 1, 2 or 3 for 2D or 3D plots
grid2Dplot = 'PhTh'; % Can be DirCos, TrueView, PhTh, etc
coorPlot = 'spherical';
output_1and2D = 'E2';

rotHandle = @rotGRASP;
rotAng = deg2rad([45,45,0]);

% Read a test field
pathName = 'Farfield Source [1]';
FF = FarField.readCSTffs(pathName);


handle2Dgrid = str2func(['grid2',grid2Dplot]);
handleCoor = str2func(['coor2',coorPlot]);
% Plot the original version
if plotDim == 3
    FF.plot('plotType','3D')
elseif plotDim == 2
    FF = handle2Dgrid(FF);
    FF = handleCoor(FF,false);
    FF.plot('plotType','2D','output',output_1and2D,'outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40), hold on
elseif plotDim == 1
    FF = handleCoor(FF,false);
    FF.plotPrincipleCuts('output',output_1and2D);
end

% Rotate the field
FFr = FF.rotate(rotHandle,rotAng);
if plotDim == 3
    figure
    FFr.plot('plotType','3D')
elseif plotDim == 2
    figure
    FFr = handle2Dgrid(FFr);
    FFr = handleCoor(FFr,false);
    FFr.plot('plotType','2D','output',output_1and2D,'outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40), hold on
elseif plotDim == 1
    FFr = handleCoor(FFr,false);
    FFr.plotPrincipleCuts('output',output_1and2D);
end

