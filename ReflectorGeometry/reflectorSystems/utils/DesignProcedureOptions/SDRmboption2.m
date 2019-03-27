function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRmboption2(Dm,F,th_e, sigma, Df)
f = SDReq28(F,Dm,Df,sigma,th_e);
Lm = SDReq29(F,f);
Ds = SDReq25(F,Df,f);
Ls = SDReq27(F,Dm,f,sigma);
a = SDReq6(Ls,f);
end