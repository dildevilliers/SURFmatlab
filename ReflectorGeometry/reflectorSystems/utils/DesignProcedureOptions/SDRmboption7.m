function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRmboption7(Dm,Ls,th_e, sigma, Df)
[F1,F2] = SDReq34(Df,Dm,sigma,Ls,th_e);

f1 = SDReq35(Ls,th_e,Dm,F1,sigma);
Ds1 = SDReq25(F1,Df,f1);
Lm1 = SDReq29(F1,f1);
a1 = SDReq6(Ls,f1);

f2 = SDReq35(Ls,th_e,Dm,F2,sigma);
Ds2 = SDReq25(F2,Df,f2);
Lm2 = SDReq29(F2,f2);
a2 = SDReq6(Ls,f2);

[Dm_r1,F_r1,Lm_r1,Ds_r1,Ls_r1,a_r1,f_r1,th_e_r1] = SDRmboption1(Dm,F1,Lm1,sigma, Df);
[Dm_r2,F_r2,Lm_r2,Ds_r2,Ls_r2,a_r2,f_r2,th_e_r2] = SDRmboption1(Dm,F2,Lm2,sigma, Df);

r1_ev = abs(Dm_r1-Dm) + abs(F_r1-F1) + abs(Lm_r1-Lm1) + abs(Ds_r1-Ds1) + abs(Ls_r1-Ls) + abs(a_r1-a1) + abs(f_r1-f1) + abs(th_e_r1-th_e);
r2_ev = abs(Dm_r2-Dm) + abs(F_r2-F2) + abs(Lm_r2-Lm2) + abs(Ds_r2-Ds2) + abs(Ls_r2-Ls) + abs(a_r2-a2) + abs(f_r2-f2) + abs(th_e_r2-th_e);

if r1_ev < r2_ev
    a = a1;
    Ds = Ds1;
    f = f1;
    F = F1;
    Lm = Lm1;
else
    a = a2;
    Ds = Ds2;
    f = f2;
    F = F2;
    Lm = Lm2;
end
end