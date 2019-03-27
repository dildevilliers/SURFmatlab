function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRmboption6(Dm,Lm,th_e, sigma, Df)
F = SDReq33(Df,Dm,sigma,Lm,th_e);
f = SDReq10(Lm,F);
Ds = SDReq25(F,Df,f);
Ls = SDReq27(F,Dm,f,sigma);
a = SDReq6(Ls,f);
end