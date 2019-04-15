function [thHat,phHat,rHat] = unitVectorsSpherical(th,ph)
% Returns the 3 spherical coordinate system unit vectors as a function of
% input angles th and ph. Output vectors in pnt3D class

tHatx = cos(th).*cos(ph);
tHaty = cos(th).*sin(ph);
tHatz = -sin(th);
thHat = pnt3D(tHatx,tHaty,tHatz);

pHatx = -sin(ph);
pHaty = +cos(ph);
pHatz = zeros(size(ph));
phHat = pnt3D(pHatx,pHaty,pHatz);

rHatx = sin(th).*cos(ph);
rHaty = sin(th).*sin(ph);
rHatz = +cos(th);
rHat = pnt3D(rHatx,rHaty,rHatz);
