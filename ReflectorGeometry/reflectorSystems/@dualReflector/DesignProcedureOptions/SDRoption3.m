function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRoption3(Dm,F,Ls,th_e)
Lm = SDReq12(Ls,F,sigma,Dm,th_e);
Ds = SDReq13(th_e,Ls,F,sigma,Dm);
a = SDReq14(th_e,Ls,F,sigma,Dm);
f = SDReq15(th_e,Ls,F,sigma,Dm);
end