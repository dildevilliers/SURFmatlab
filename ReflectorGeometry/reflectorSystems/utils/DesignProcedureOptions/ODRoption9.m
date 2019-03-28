function [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption9(Dm, th_0, th_e, Dsx, beta, sigma)
e = ODReq5(sigma,beta,th_0);
alpha = ODReq6(e,beta);
th_U = ODReq3(e,alpha,sigma,th_e,beta);
F = ODReq26(Dm,th_U,th_0);
h = ODReq23(F,th_0);
th_L = ODReq4(h,Dm,F);
a = ODReq18(sigma,Dsx,e,beta,th_U,th_L);
f = ODReq15(a,e);
Ls = ODReq8(a,e,beta,th_0);
Lm = ODReq9(a,e,beta,th_0,h);
end