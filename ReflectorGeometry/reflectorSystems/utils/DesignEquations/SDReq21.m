function f = SDReq21(th_e,F,sigma,Dm,Ds)
Y1 = tan(th_e./2).*F;
Y2 = 8.*F.*Dm + sigma.*sin(th_e).*(16.*F.^2 + Dm.^2);
f = (Ds.*Y2.*(sigma.*Dm - 4.*Y1))./(32.*sin(th_e).*F.*(sigma.*Dm + 4.*Y1).*Dm);
end