function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRoption7(Dm,Ds,Ls,th_e)
F = SDReq24(th_e,Ls,Ds,Dm,sigma);
f = SDReq15(th_e,Ls,F,sigma,Dm);
a = SDReq14(th_e,Ls,F,sigma,Dm);
Lm = SDReq12(Ls,F,sigma,Dm,th_e);
end