% CBFP Test Script

% Read in Farfields
% Options:  1) [1 x 1]:      FarField object with multiple frequencies - FFobj1
%           2) [1 x nFFobj]: Multiple FarField objects defined at multiple frequencies (Frequency samples must be the same for each FarField object) - FFobj2
%           3) [1 x nFFobj]: Multiple FarField objects with a single frequency sample. Frequency treated as a normal parameter. - FFobj3
%

FFobj1 = FarField.readCSTffs('C:\Users\21502382\Documents\GitHub\SURFmatlab\testScripts\Farfield Source [1]');

FFobj = FFobj1;

%% DO CBFP DECOMPOSITION
% Options: 1) Use all FarFields in the FarField object to construct the CBFP object - CBFPobj1
%          2) Specify tolerance, basis function indices or total number of basis functions to use (selection from 1 to nBasis) - CBFPobj2

% Selection
tol1 = 1e-100; % Specify largest singular value to keep. Default = 1e-100
iFM1 = []; % Ignore indices selection
nBasis1 = 9; % Use the first 9 out of 11 FarField patterns for the CBFP decomposition

% Call Function
CBFPobj1 = CBFP(FFobj);
CBFPobj2 = CBFP(FFobj,tol1,iFM1,nBasis1);

%% FIELD CONSTRUCTION FROM BASIS FUNCTIONS
% Options: 1) Re-construct FarFields used to generate the CBFP basis fuctions in the CBFP object - FF1
%          2) Re-construct FarFields only at the weights specified in the coefficient matrix [W] 
%               [W] options: 1) [nBasis x nSupportPoints]
%                            2) [nBasis x nSupportPoints x nFrequencies]
%
%          3) Select which basis functions to re-construct the FarFields with

% Selection
tol2 = 1e-100; % Specify largest singular value to keep. Default = 1e-100
iBasis2 = []; % Ignore indices selection
nBasis2 = 9; % Use the first 9 out of 11 FarField patterns for the CBFP decomposition

% Call Function
FF1 = CBFP.expansion2FarField(CBFPobj1);

% FF2 = expansion2FarField(CBFPobj,W);
% 
% FF3 = expansion2FarField(CBFPobj,W,tol2,iBasis2,nBasis2);

%% Plots
% Plot first (most significant) basis function
% Plot first FarField pattern

% BASIS FUNCTION
basis1 = CBFPobj1.basis(1,1);

% % 1D plot
% basis1.plot('dynamicRange_dB',40,'plotType','cartesian','output','Directivity','step',1);

% % 2D plot
% figure('Name','Basis Function 2D - E1')
% basis1.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40);
% 
% 3D plot
figure('Name','Basis Function_1')
basis1.plot;
% 
% % RECONSTRUCTED FARFIELD
% % 1D plot
% FF1.plotPrincipleCuts('dynamicRange_dB',40,'plotType','cartesian','output','Directivity');
% 
% % 2D plot
% figure('Name','FarField 2D - E1')
% FF1.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40);
% 
% 3D plot
figure('Name','FarField')
FF1.plot;

%% Interpolate and Plot weights


