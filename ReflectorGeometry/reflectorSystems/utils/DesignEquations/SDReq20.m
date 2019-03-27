function a = SDReq20(th_e,F,sigma,Dm,Ds)
Y2 = 8.*F.*Dm + sigma.*sin(th_e).*(16.*F.^2 + Dm.^2);
a = (Ds.*Y2)./(32.*sin(th_e).*Dm.*F);
end