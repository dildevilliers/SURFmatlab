close all
clear all

exNumber = 1;

th_ext = deg2rad(20);
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
        DR = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext);
    case 2
        % Granet Dual Reflector Ex2
        Dm = 100;
        Lm = 109.249;
        th_e = deg2rad(11.9131);
        Ls = 41.2498;
        th_0 = deg2rad(-39.0356);
        beta = deg2rad(5.4);
        sigma = 1;
        DR = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext);
    case 3
        % Granet Dual Reflector Ex3
        Dm = 45;
        Lm = 40.32365;
        th_e = deg2rad(10.32476);
        Ls = 21.04870;
        th_0 = deg2rad(-55.51708);
        beta = deg2rad(6.0);
        sigma = -1;
        DR = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext);
    case 4
         % Granet Dual Reflector Ex4
        Dm = 24;
        Lm = 34.03933;
        th_e = deg2rad(11.50497);
        Ls = 30.54596;
        th_0 = deg2rad(-53.1301);
        beta = deg2rad(5.6);
        sigma = 1;
        DR = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext);
    case 5
        Dm = 10;
        Lm = 5;
        th_e = deg2rad(25);
        Ls = 2;
        th_0 = 0;
        beta = 0;
        sigma = 1;
        DR = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext);
    case 6 
        % Granet Symmetric Cassegrain no blockage (Fig 5)
        Dm = 22;
        Lm = 1.4542;
        Ds = 2.75;
        Ls = 5.0868;
        a = 2.1244;
        f = 2.9623;
        th_e = deg2rad(14);
        th_0 = 0;
        beta = 0;
        sigma = -1;
        DR = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext);
        if abs(DR.a - a) > 1e-4 || abs(DR.f - f) > 1e-4 
            warning('Case 6 broken')
        end
    case 7
        % Granet Symmetric Gregorian no blockage (Fig 6)
        Dm = 10;
        F = 3.0928;
        Lm = 1;
        Ds = 1.1898;
        Ls = 2.5;
        a = 1.4535;
        f = 1.0464;
        th_e = deg2rad(15);
        th_0 = 0;
        beta = 0;
        sigma = 1;
        DR = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext);
        if abs(DR.a - a) > 1e-4 || abs(DR.f - f) > 1e-4  
            warning('Case 7 broken')
        end
    case 8
        % Granet Symmetric Cassegrain minimum blockage (Fig 7)
        Dm = 10;
        F = 5;
        Lm = 1;
        Ds = 1.25;
        Ls = 3.4024;
        a = 1.4024;
        f = 2;
        th_e = deg2rad(10.0369);
        th_0 = 0;
        beta = 0;
        sigma = -1;
        DR = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext);
        if abs(DR.a - a) > 1e-4 || abs(DR.f - f) > 1e-4 || abs(DR.F - F) > 1e-4 
            warning('Case 8 broken')
        end
end
figure
DR.plot3D(10000,[1,1,1])
figure
DR.plot

% Px1 = DR.SR.coor.origin.pointMatrix + DR.SR.coor.x_axis.*3;
% Px2 = DR.SR.coor.origin.pointMatrix + DR.SR.coor.x_axis.*0;
% Pz1 = DR.SR.coor.origin.pointMatrix + DR.SR.coor.z_axis.*3;
% Pz2 = DR.SR.coor.origin.pointMatrix + DR.SR.coor.z_axis.*0;
% plot([Px1(1),Px2(1)],[Px1(3),Px2(3)],'k')
% plot([Pz1(1),Pz2(1)],[Pz1(3),Pz2(3)],'r')

%%
figure
DR.plotRayTrace
pathLengthStruct = DR.getPathLength
[phPRmask,thPRmask] = meshgrid(linspace(0,2*pi,101),linspace(0,pi,101));
[FFM,MaskPointing,M] = DR.getMask([phPRmask(:),thPRmask(:)],'SR');

figure
FFM = FFM.grid2TrueView;
plot(FFM,'output','Directivity','outputType','mag','plotType','2D','scaleMag','lin','step',1,'norm',1)
% Also plot the mask pointing
[x,y,z] = deal(zeros(size(MaskPointing)));
for mm = 1:1:length(MaskPointing)
    endPoints = MaskPointing(mm).origin.addVect(MaskPointing(mm).z_axis);
    x(mm) = endPoints.x;
    y(mm) = endPoints.y;
    z(mm) = endPoints.z;
end
figure
plot3(x(M==2),y(M==2),z(M==2),'k*'), grid on, hold on
plot3(x(M==1),y(M==1),z(M==1),'b.')
plot3(x(M==0),y(M==0),z(M==0),'r.')
axis equal







