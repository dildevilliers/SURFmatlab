% Test the symmetry implementation of the FarField class
close all
clear all

%% Get something to work with
pathName = 'Farfield Source [1]';
FF = FarField.readCSTffs(pathName);
FF = FF.coor2Ludwig3(false);
FF = FF.setXrange('sym');
FF = FF.currentForm2Base;

% Plot it
scaleMag = 'lin';
outputType = 'real';
output = 'E2';
plotGridHandle = @grid2PhTh;
FFplot = plotGridHandle(FF);
FFplot.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)

%% Get a version symmetric about x
% phXind = FF.ph >= pi;
% phXind = FF.ph >= -eps;
% FFx = FarField(FF.ph(phXind),FF.th(phXind),FF.E1(phXind,:),FF.E2(phXind,:),FF.E3(phXind,:),FF.freq,FF.Prad./2,FF.radEff,FF.coorSys,FF.polType,'PhTh',FF.freqUnit);
% FFx = FFx.setSymmetryXZ('electric');
% FFx = FFx.setXrange('sym');
% FFx = plotGridHandle(FFx);
% figure
% FFx.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
% % Fix to full range and plot
% FFxFull = FFx.getFullPattern;
% FFxFull = plotGridHandle(FFxFull);
% figure
% FFxFull.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
% 
% FFdelta = FFxFull - FF;
% figure
% FFdelta.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType','mag','scaleMag',scaleMag)
% [eEx] = FFdelta.norm

%% Get a version symmetric about y
phYind = FF.ph >= -pi/2-eps & FF.ph <= pi/2+eps;
FFy = FarField(FF.ph(phYind),FF.th(phYind),FF.E1(phYind,:),FF.E2(phYind,:),FF.E3(phYind,:),FF.freq,FF.Prad./2,FF.radEff,FF.coorSys,FF.polType,'PhTh',FF.freqUnit);
FFy = FFy.setSymmetryYZ('magnetic');
FFy = FFy.setXrange('sym');
FFy = plotGridHandle(FFy);
figure
FFy.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
% Fix to full range and plot
FFyFull = FFy.getFullPattern;
FFyFull = plotGridHandle(FFyFull);
figure
FFyFull.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)

FFdelta = FFyFull - FF;
figure
FFdelta.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType','mag','scaleMag',scaleMag)
[eEy] = FFdelta.norm
