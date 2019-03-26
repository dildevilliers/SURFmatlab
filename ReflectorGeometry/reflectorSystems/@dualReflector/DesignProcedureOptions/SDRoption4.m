function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRoption4(F,Ds,Ls,th_e)
Dm = SDReq16(th_e,Ls,F,sigma,Ds);
Lm = SDReq12(Ls,F,sigma,Dm,th_e);
f = SDReq15(th_e,Ls,F,sigma,Dm);
a = SDReq14(th_e,Ls,F,sigma,Dm);
end