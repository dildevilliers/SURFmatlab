function Ls = SDReq27(F,Dm,f,sigma,th_e)
Ls = (2.*sigma.*Dm.*f)./(sigma.*Dm - 4.*F.*tan(th_e./2));
end