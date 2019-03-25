function [F,h,th_U,th_L,e,a,f,Dsx,Dsy,alpha,Lm,dSRPR,dFPR,Lt,Ht,C_SR] = ODRoption8(Dm, th_0, th_e, Ls, beta)
e = ODReq5(sigma,beta,th_0);
alpha = ODReq6(e,beta);
a = ODReq22(Ls,e,beta,th_0);
f = ODReq15(a,e);
th_U = ODReq3(e,alpha,sigma,th_e,beta);
F = ODReq26(Dm,th_U,th_0);
h = ODReq23(F,th_0);
th_L = ODReq4(h,Dm,F);
Dsx = ODReq27(a,sigma,e,beta,th_U,th_L);
Lm = ODReq9(a,e,beta,th_0,h);
dSRPR = ODReq11(h,Dm,a,sigma,e,th_U,beta,th_L);
dFPR = ODReq10(h,Dm,f,beta);
Lt = ODReq12(a,sigma,e,th_L,beta,th_U,Dm,F,h);
Ht = ODReq13(h,Dm,a,sigma,e,th_L,beta,th_U);
C_SR = ODReq38(a,sigma,e,beta,th_U,th_L,th_e,alpha,f);
Dsy = ODReq39(a,e,th_e,alpha);
end