function C_SR = ODReq38length(F0P1,F0P2,a,sigma,th_e,alpha,f)
Cx = (F0P1.*sin(alpha + sigma.*th_e) + F0P2.*sin(alpha - sigma.*th_e))./2;
Cy = 0;
Cz = a.*sqrt(1 + Cx.^2/(f.^2 - a.^2)) - f;
C_SR = Pnt3D(Cx,Cy,Cz);
end