function [Dm, th_0, th_e, Ls, beta] = ODRoption2(Dm, F, h, Ls, beta)
th_0 = ODReq1(h,F);
th_U = ODReq2(h,Dm,F);
th_L = ODReq4(h,Dm,F);
e = ODReq5(sigma,beta,th_0);
alpha = ODReq6(e,beta);
th_e = ODReq7(sigma,e,th_U,beta,alpha);
end