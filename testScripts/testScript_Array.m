% Test script for the suite of Array classes

close all
clear all

%%
fRF = 1.57542e9;
th_in = deg2rad(0);
ph_in = deg2rad(0);
Ps = -50; % in dBm
phaseSig = deg2rad(0);

Nant = 4;      % Number of elements
d_lam = 2;     % Element spacing in wavelengths 
d = d_lam*(physconst('lightspeed')/fRF);    % Element spacing in m
x = (-(Nant-1)/2:(Nant-1)/2).*d;       % Uniform linear array along x-axis
r = pnt3D(x,0,0);
ampErrors = [1,1,1,1];
phaseErrors_deg = [0,0,0,0];
channelPhasors = ampErrors.*exp(1i.*deg2rad(phaseErrors_deg));

Pn = -110; % in dBm
LNAGain = 20;
IFGain = 60;
fIF = 4.096e6;
fLO = fRF - fIF; 

fSamp = 4*fIF;          % Sample rate in Hz
Nt = 8000;                % Number of time samples
delT = 1/fSamp;
t0 = 0;
t = t0:delT:(t0+delT*(Nt-1));

Nbits = 1;
A = ArrayADC(Nbits);

S1 = PlaneWaveSignal('compExp',fRF,th_in,ph_in,Ps,phaseSig);
S2 = PlaneWaveSignal('compExp',fRF,th_in + deg2rad(40),ph_in,Ps,phaseSig+pi/2);
S3 = PlaneWaveSignal('noise',fRF,th_in + deg2rad(-20),ph_in,Ps,phaseSig);
% S = [S1,S2,S3];
S = S1;
% E = ArrayElements(r);
% R = ArrayReceiver(Pn,LNAGain,IFGain,fLO);
% portSigMat = E.portSignals(S,t);
% sn = R.sigRec(portSigMat,t);
% % b = A.getBinarySignal(real(sn));
% % x = A.getDigitalSignal(b);
% x = A.ADC(real(sn));
% plot(t,x(1,:))
% 
% plot(t,real(sn(1,:)),'k'), grid on, hold on

arraySys = ArraySystem(r,channelPhasors,Pn,LNAGain,IFGain,fLO,fSamp,Nt,Nbits);
arraySys.plotPortSignal(S);
x = arraySys.getPortSignal(S);

arrayDBE = ArrayDBE(arraySys);
figure
arrayDBE.plotPortPSD(x,1,4)

% Test the scanning
Nscan = 501;
th_scan_deg = linspace(-60,60,Nscan);
% ph_scan_deg = ones(size(th_scan_deg))*0;    % Use ph = 0 for array along x-axis
ph_scan_deg = 0;    % Use ph = 0 for array along x-axis
% y = arrayDBE.scanBeam(fRF,deg2rad(th_scan_deg),deg2rad(ph_scan_deg),x);
figure
arrayDBE.plotScanBeam(fRF,deg2rad(th_scan_deg),deg2rad(ph_scan_deg),x);
