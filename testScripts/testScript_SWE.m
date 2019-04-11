FF1 = FarField.readFEKOffe([pwd,'\patch']);
FF2 = FarField.readCSTffs('Farfield Source [1]');

thlims = [0 pi/2];
phlims = [0 2*pi];

%Create SWE object with target FarField and MMAX/NMAX numbers
swe1 = SWE(FF1,[4 6],'MNmax');

%Create SWE object with target FarField and minimum sphere radius
swe2 = SWE(FF1,0.01,'r0');

%Create SWE object with set of angle vectors and MMAX/NMAX numbers
th = FF1.y; %linspace(0,pi,37);
ph = FF1.x; %linspace(0,2*pi,37);
%Select angular sector
isec = intersect( intersect( find( th >= thlims(1)),find( th <= thlims(2))) , intersect( find( ph >= phlims(1)),find( ph <= phlims(2))) );%% && th <= thlims(2) && ph >= phlims(1) && ph <= phlims(2) );

thsec = th(isec);
phsec = ph(isec);
swe3 = SWE([inf.*ones(size(thsec)) thsec phsec], [10 10], 'MNmax');

% %Create SWE object with target FarField (multiple frequencies) and MMAX/NMAX numbers
% swe4 = SWE(FF2,[6 4],'MNmax');
% %Create SWE object with target FarField (multiple frequencies) and MMAX/NMAX numbers
% swe5 = SWE(FF2,[10 10],'MNmax');
% %Create SWE object with target FarField (multiple frequencies), MMAX/NMAX numbers and an explicit frequency setting
% swe6 = SWE(FF2,[6 4],'MNmax',1e9);
% 
% %Perform SWE expansion on patch farfield
% [Qj1,P1] = SWE.farField2Expansion(swe1,FF1);
% [Qj2,P2] = SWE.farField2Expansion(swe2,FF1);
[Qj3,P3] = SWE.farField2Expansion(swe3,FF1);
[Qj4,P4] = SWE.farField2Expansion(swe4,FF2);
[Qj5,P5] = SWE.farField2Expansion(swe5,FF2);
[Qj6,P6] = SWE.farField2Expansion(swe6,FF2);

%Reconstruct patch farfield from SWE expansion coefficients
FF_swe1 = SWE.expansion2FarField(swe1,Qj1,P1);
FF_swe2 = SWE.expansion2FarField(swe2,Qj2,P2);
FF_swe3 = SWE.expansion2FarField(swe3,Qj3,P3);
FF_swe4 = SWE.expansion2FarField(swe4,Qj4,P4);
FF_swe5 = SWE.expansion2FarField(swe5,Qj5,P5);
FF_swe6 = SWE.expansion2FarField(swe6,Qj6,P6);

%Plot field differences between target and SWE-reconstructed fields
%FEKO farfield (1 frequency)
figure
subplot(3,1,1)
plot(FF1-FF_swe1{1},'plotType','2D','outputType','mag','step',1,'showGrid',true);
subplot(3,1,2)
plot(FF1-FF_swe2{1},'plotType','2D','outputType','mag','step',1,'showGrid',true);
subplot(3,1,3)
plot(FF1-FF_swe3{1},'plotType','2D','outputType','mag','step',1,'showGrid',true);

%CST farfield (11 frequencies)
figure
subplot(3,1,1)
plot(FF2-FF_swe4{1},'plotType','2D','outputType','mag','step',1,'showGrid',true);
subplot(3,1,2)
plot(FF2-FF_swe5{1},'plotType','2D','outputType','mag','step',1,'showGrid',true);
subplot(3,1,3)
plot(FF2-FF_swe6{1},'plotType','2D','outputType','mag','step',1,'showGrid',true);

