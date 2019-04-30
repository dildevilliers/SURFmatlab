clear all
close all

exNumber = 2;

th_ext = deg2rad(20); % Extension angle in (rad)
symFact_ext = 1;

switch exNumber
    case 0
        % The default - MeerKAT with 20 degree extension
        DR = dualReflector;
    case 1
        % Granet Dual Reflector Ex1
        Dm = 100;
        Lm = 107.772;
        th_e = deg2rad(11.8767);
        Ls = 28.0096;
        th_0 = deg2rad(-40.608);
        beta = deg2rad(10.1);
        sigma = -1;
    case 2
        % Granet Dual Reflector Ex2
        Dm = 100;
        Lm = 109.249;
        th_e = deg2rad(11.9131);
        Ls = 41.2498;
        th_0 = deg2rad(-39.0356);
        beta = deg2rad(5.4);
        sigma = 1;
    case 3
        % Granet Dual Reflector Ex3
        Dm = 45;
        Lm = 40.32365;
        th_e = deg2rad(10.32476);
        Ls = 21.04870;
        th_0 = deg2rad(-55.51708);
        beta = deg2rad(6.0);
        sigma = -1;
    case 4
         % Granet Dual Reflector Ex4
        Dm = 24;
        Lm = 34.03933;
        th_e = deg2rad(11.50497);
        Ls = 30.54596;
        th_0 = deg2rad(-53.1301);
        beta = deg2rad(5.6);
        sigma = 1;
    case 5
        Dm = 10;
        Lm = 5;
        th_e = deg2rad(25);
        Ls = 2;
        th_0 = 0;
        beta = 0;
        sigma = 1;
end

DR = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext);
%DR = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma);
%exportDualReflectorToTOR(obj,fullpathName,freqValGhz,prefixName)
DR.exportDualReflectorToTOR('C:\testing_area\fullDoubleExampleTor.tor',6,'dubble');
%GRASP_angle = (2.*DR.f)/DR.F
figure
DR.plot3D(10000,[1,1,1])
figure
DR.plot
figure
DR.plotRayTrace
