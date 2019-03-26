function Ls = SDReq27(F,Dm,f,sigma)
Ls = (2.*sigma.*Dm.*f)./(sigma.*Dm - 4.*F.*tan(th_e./2));
end