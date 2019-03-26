function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRoption1(Dm,Lm,Ls,th_e)
f = SDReq4(Ls,sigma,Dm,Lm,th_e);
F = SDReq5(Lm,f);
a = SDReq6(Ls,f);
Ds = SDReq7(Ls,f,th_e,sigma,F,Dm);
end