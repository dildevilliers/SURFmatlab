function [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption7(Dm, th_0, dFPR, Ls, beta, sigma)
e = ODReq5(sigma,beta,th_0);
alpha = ODReq6(e,beta);
a = ODReq22(Ls,e,beta,th_0);
f = ODReq15(a,e);
h = ODReq24(dFPR,Dm,f,beta);
F = ODReq25(h,th_0);
th_U = ODReq2(h,Dm,F);
th_L = ODReq4(h,Dm,F);
th_e = ODReq7(sigma,e,th_U,beta,alpha);
Lm = ODReq9(a,e,beta,th_0,h);
end