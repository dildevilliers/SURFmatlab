
%create farfield object
FF = FarField.readFEKOffe([pwd,'\infdipole']);


%2D plot of mag(E1), 1 degree step, 40dB dynamic range and visible grid 
figure
FF.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40);
%2D plot of mag(E2), 1 degree step, 40dB dynamic range and visible grid 
figure
FF.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40);

%1D Cartesian plot of directivity
figure
FF.plot('plotType','cartesian','output','Directivity','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40,'cutConstant','ph');
