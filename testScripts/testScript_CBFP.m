% testScript_CBFP

%% Create FarField Objects
% Options:  1) [1x1]: FarField object with a single freq - FFobj1
%           2) [1x1]: FarField object with multiple freqs - FFobj2
%           3) [1xnFFobjs]: An array of FarField objects with a single freq. Frequency treated as a normal parameter. - FFobj1s
%           4) [1xnFFobjs]: An array of FarField objects with multiple freqs. (Frequency samples must be the same for each FarField object) - FFobj2s

FFobj1 = FarField();
FFobj2 = FarField.readCSTffs('Farfield Source [1]');

FFobj1s = [FFobj1,FFobj1,FFobj1,FFobj1];
FFobj2s = [FFobj2,FFobj2,FFobj2,FFobj2,FFobj2];

%% DO CBFP DECOMPOSITION
% CBFPobj = CBFP(FFobj,tol,iFM,iBasis)
% 3 Cases:  1) Expansion across freq (CBFPobj1)
%           2) Freq treated as a normal parameter (CBFPobj2)
%           3) Geometric variation per freq (CBFPobj3)

CBFPobj1 = CBFP(FFobj2);
CBFPobj2 = CBFP(FFobj1s);
CBFPobj3 = CBFP(FFobj2s);

% SELECTION
tol = 1e-10;
iFM1.f = [1,3,5,7,9,11];
iFM2.x = [2,4];
iBasis = [1,2,3,4];

CBFPobj_a = CBFP(FFobj1s,tol); % For this particular example: Truncates basis functions from 4 to 1
CBFPobj_b = CBFP(FFobj2,[],iFM1); % Which FarFields to use across freq
CBFPobj_c = CBFP(FFobj2s,[],iFM2); % Which FarFields to use across geometric variation
CBFPobj_d = CBFP(FFobj2,[],[],iBasis); % How many basis functions to return

%% GET COEFFICIENTS FROM A TARGET FIELD
% Wout = farField2Coeffs(FFobj,CBFPobj,tol,iBasis)

Fobj = getFi(FFobj2,7); % get 7nth freq from FFobj2
Fobjs(1,1) = getFi(FFobj2,3);
Fobjs(1,2) = getFi(FFobj2,4); % [1x2] array of FarFields

W1 = CBFP.farField2Coeffs(Fobj,CBFPobj1); % 1 Target Field as an input
W2 = CBFP.farField2Coeffs(Fobjs,CBFPobj1); % Multiple Target Fields as an array input 
W3 = CBFP.farField2Coeffs(FFobj1,CBFPobj2);
W4 = CBFP.farField2Coeffs(Fobj,CBFPobj3); 

% SELECTION
iBasis = (1:3);
W5 = CBFP.farField2Coeffs(Fobj,CBFPobj1,[],iBasis); % Get coeffs of the first 3 basis functions

%% FIELD CONSTRUCTION FROM BASIS FUNCTIONS
% FFobj = coeffs2FarField(CBFPobj,W,tol,iBasis,freqIndex)

FF1 = CBFP.coeffs2FarField(CBFPobj1); % Reconstruct FarFields used to generate basis functions
FF2 = CBFP.coeffs2FarField(CBFPobj1,W1); 

% SELECTION
FF3 = CBFP.coeffs2FarField(CBFPobj2,[],[],iBasis); % Reconstruct FarFields using basis functions 1,2 & 3
W6 = {CBFPobj3.coeffs{1,1}(:,7)}; 
FF4 = CBFP.coeffs2FarField(CBFPobj3,W6,[],[],7); % CBFP Decomposition Case 3: At freq slice 7, reconstruct FarField using geomtrically varying basis funcs

%% Plots

% FarField obj Options
% FF1 ...  FF4
% FFobj1 ... FFobjs(1,b)
% CBFPobj.basis{a,b}

basisFn1 = CBFPobj1.basis{1,1};

plotObj = basisFn1;

% 1D plot
plotObj.plot('output','Directivity','outputType','mag','plotType','cartesian','scaleMag','dB','norm',1,...
    'step',1,'dynamicRange_dB',50,'freqIndex',1,'cutValue',0)

% 2D plot
figure()
plotObj.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40);

% 3D plot
figure()
plotObj.plot;
