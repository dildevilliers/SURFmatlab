function Ds = SDReq11(th_e,Lm,F,sigma,Dm)
X1 = -16.*sin(th_e).*Dm.*F.*(Lm - F).*(sigma.*Dm + 4.*F.*tan(th_e./2));
X2 = 8.*F.*Dm.*(sigma.*Dm - 4.*F.*tan(th_e/2));
X3 = (Dm.^2 + 16.*F.^2).*sin(th_e).*(Dm - 4.*sigma.*F.*tan(th_e./2));
Ds = (X1)./(X2 + X3);
end