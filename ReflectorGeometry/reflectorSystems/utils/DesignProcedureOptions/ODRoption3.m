function [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption3(Dm, F, h, DFPR, beta, sigma)
th_0 = ODReq1(h,F);
th_U = ODReq2(h,Dm,F);
th_L = ODReq4(h,Dm,F);
e = ODReq5(sigma,beta,th_0);
alpha = ODReq6(e,beta);
th_e = ODReq7(sigma,e,th_U,beta,alpha);
f = ODReq16(dFPR,h,Dm,beta);
a = ODReq17(f,e);
Dsx = ODReq27(sigma,a,e,beta,th_U,th_L);
Lm = ODReq9(a,e,beta,th_0,h);
Ls = ODReq8(a,e,beta,th_0);
end