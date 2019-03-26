function [Dm, th_0, th_e, Ls, beta] = ODRoption1(Dm, F, h, Dsx, beta)
th_0 = ODReq1(h,F);
th_U = ODReq2(h,Dm,F);
th_L = ODReq4(h,Dm,F);
e = ODReq5(sigma,beta,th_0);
alpha = ODReq6(e,beta);
th_e = ODReq7(sigma,e,th_U,beta,alpha);
a = ODReq18(sigma,Dsx,e,beta,th_U,th_L);
f = ODReq15(a,e);
Ls = ODReq8(a,e,beta,th_0);
end