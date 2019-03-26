function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRoption5(Lm,Ds,Ls,th_e)
f = SDReq17(th_e,Ls,Ds);
Dm = SDReq18(th_e,Lm,sigma,Ls,f);
F = SDReq5(Lm,f);
a = SDReq6(Ls,f);
end