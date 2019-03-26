function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRmboption1(Dm,F,Lm,sigma)
Df = 1;
f = SDReq10(Lm,F);
Ds = SDReq25(F,Df,f);
th_e = SDReq26(F,Dm,Ds,f,sigma);
Ls = SDReq27(F,Dm,f,sigma,th_e);
a = SDReq6(Ls,f);
end