function C_SR = ODReq38(a,sigma,e,beta,th_U,th_L,th_e,alpha,f)
F0P1 = ODReq30_F0P1(a,sigma,e,beta,th_L);
F0P2 = ODReq30_F0P2(a,sigma,e,beta,th_U);

Cx = (F0P1.*sin(alpha + sigma.*th_e) + F0P2.*sin(alpha - sigma.*th_e))./2;
Cy = 0;
Cz = a.*sqrt(1 + Cx.^2/(f.^2 - a.^2)) - f;
C_SR = pnt3D(Cx,Cy,Cz);
end