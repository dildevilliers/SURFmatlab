function a = SDReq8(Lm,F,sigma,Dm,th_e)
a = (-1/2).*( ((Lm-F).*(sigma.*Dm + 4.*F.*tan(th_e./2)))./(sigma.*Dm - 4.*F.*tan(th_e./2)) );
end