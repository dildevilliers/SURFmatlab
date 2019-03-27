function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRmboption7(Dm,Ls,th_e, sigma, Df)
F = SDReq34(Df,Dm,sigma,Ls,th_e);
f = SDReq35(Ls,th_e,Dm,F,sigma);
Ds = SDReq25(F,Df,f);
Lm = SDReq29(F,f);
a = SDReq6(Ls,f);
end