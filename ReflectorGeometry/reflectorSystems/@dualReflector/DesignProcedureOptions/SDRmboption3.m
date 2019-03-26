function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRmboption3(Dm,F,Ds)
f = SDReq30(F,Df,Ds);
Lm = SDReq29(F,f);
th_e = SDReq26(F,Dm,Ds,f,sigma);
Ls = SDReq27(F,Dm,f,sigma);
a = SDReq6(Ls,f);
end