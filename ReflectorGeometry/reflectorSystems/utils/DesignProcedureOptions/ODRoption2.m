function [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption2(Dm, F, h, Ls, beta, sigma)
th_0 = ODReq1(h,F);
th_U = ODReq2(h,Dm,F);
th_L = ODReq4(h,Dm,F);
e = ODReq5(sigma,beta,th_0);
alpha = ODReq6(e,beta);
th_e = ODReq7(sigma,e,th_U,beta,alpha);
a = ODReq22(Ls,e,beta,th_0);
f = ODReq15(a,e);
Dsx = ODReq27(a,sigma,e,beta,th_U,th_L);
Lm = ODReq9(a,e,beta,th_0,h);
end