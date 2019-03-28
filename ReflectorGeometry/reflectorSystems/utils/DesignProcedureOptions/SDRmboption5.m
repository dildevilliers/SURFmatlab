function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRmboption5(Dm,Ds,th_e, sigma, Df)
F = SDReq32(Df,Dm,sigma,Ds,th_e);
f = SDReq30(F,Df,Ds);
Lm = SDReq29(F,f);
Ls = SDReq27(F,Dm,f,sigma);
a = SDReq6(Ls,f);
end