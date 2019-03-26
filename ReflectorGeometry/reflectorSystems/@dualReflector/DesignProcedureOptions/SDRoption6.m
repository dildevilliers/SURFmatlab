function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRoption6(Dm,F,Ds,th_e)
Ls = SDReq19(th_e,F,sigma,Dm,Ds);
a = SDReq20(th_e,F,sigma,Dm,Ds);
f = SDReq21(th_e,F,sigma,Dm,Ds);
Lm = SDReq22(th_e,F,sigma,Dm,Ds);
end