function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRoption2(Dm,F,Lm,th_e,sigma)
a = SDReq8(Lm,F,sigma,Dm,th_e);
Ls = SDReq9(Lm,F,sigma,Dm,th_e);
f = SDReq10(Lm,F);
Ds = SDReq11(th_e,Lm,F,sigma,Dm);
end