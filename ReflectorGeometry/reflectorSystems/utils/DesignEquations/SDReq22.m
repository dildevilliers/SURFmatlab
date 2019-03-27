function Lm = SDReq22(th_e,F,sigma,Dm,Ds)
Y1 = tan(th_e./2).*F;
Y3 = 8.*sigma.*Dm.^2.*F.*(Ds - 2.*sin(th_e).*F) - 32.*tan(th_e./2).*F.^2.*Dm.*(Ds + 2.*sin(th_e).*F) + Ds.*sin(th_e).*(16.*F.^2 + Dm.^2).*(Dm - 4.*sigma.*tan(th_e./2).*F);
Lm = (-1/16).*((Y3)./(Dm.*sin(th_e).*F.*(sigma.*Dm + 4.*Y1)));
end