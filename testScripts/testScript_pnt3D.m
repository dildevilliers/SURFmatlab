% Test script for the pnt3D class

close all
clear all

% M = 1000;
% N = 1000;
% x = ones(M,N).*3;
% y = ones(M,N).*2;
% z = ones(M,N).*-3;
x = [1:1:5].';
y = [1:1:5].' - 4;
z = 0;

P = pnt3D(x,y,z);
P1 = translate(P,[sqrt(2),1,sqrt(2)].');
B = isequal(P1,P);
D = distanceCart(P,P1);

P2 = P1-P;

V = P2.pointMatrix;

P.plot('marker','o','markerEdgeColor',[1,0.1,1],'markerFaceColor','k','markerSize',5), hold on
P1.plot
plotLines(P1,P,'lineColor','b')
% plotVect(P,V)
