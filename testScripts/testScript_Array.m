% Test script for the suite of Array classes

close all
clear all

%%
fRF = 1.57542e9;
th_in = deg2rad(0);
ph_in = deg2rad(0);
Ps = -30; % in dBm
phaseSig = deg2rad(0);

Nant = 6;      % Number of elements
d_lam = 2;     % Element spacing in wavelengths 
d = d_lam*(physconst('lightspeed')/fRF);    % Element spacing in m
x = (-(Nant-1)/2:(Nant-1)/2).*d;       % Uniform linear array along x-axis
r = pnt3D(x,0,0);

Pn = -110; % in dBm
LNAGain = 20;
IFGain = 20;
fIF = 4.096e6;
fLO = fRF - fIF; 

fSamp = 4*fIF;          % Sample rate in Hz
Nt = 1600;                % Number of time samples
delT = 1/fSamp;
t0 = 0;
t = t0:delT:(t0+delT*(Nt-1));

S = PlaneWaveSignal('compExp',fRF,th_in,ph_in,Ps,phaseSig);
E = ArrayElements(r);
R = ArrayReceiver(Pn,LNAGain,IFGain,fLO);
portSigMat = E.portSignals(S,t);
sn = R.sigRec(portSigMat,t);

