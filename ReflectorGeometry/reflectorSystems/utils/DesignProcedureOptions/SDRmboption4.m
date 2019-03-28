function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRmboption4(Dm,Lm,Ds, sigma, Df)
f = SDReq31(Lm,Df,Ds);
F = SDReq5(Lm,f);
th_e = SDReq26(F,Dm,Ds,f,sigma);
Ls = SDReq27(F,Dm,f,sigma);
a = SDReq6(Ls,f);
end