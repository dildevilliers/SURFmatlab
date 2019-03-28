function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRmboption5(Dm,Ds,th_e, sigma, Df)
[F1,F2] = SDReq32(Df,Dm,sigma,Ds,th_e);

f1 = SDReq30(F1,Df,Ds);
Lm1 = SDReq29(F1,f1);
Ls1 = SDReq27(F1,Dm,f1,sigma);
a1 = SDReq6(Ls1,f1);

f2 = SDReq30(F2,Df,Ds);
Lm2 = SDReq29(F2,f2);
Ls2 = SDReq27(F2,Dm,f2,sigma);
a2 = SDReq6(Ls2,f2);

[Dm_r1,F_r1,Lm_r1,Ds_r1,Ls_r1,a_r1,f_r1,th_e_r1] = SDRmboption1(Dm,F1,Lm1,sigma, Df);
[Dm_r2,F_r2,Lm_r2,Ds_r2,Ls_r2,a_r2,f_r2,th_e_r2] = SDRmboption1(Dm,F2,Lm2,sigma, Df);

r1_ev = abs(Dm_r1-Dm) + abs(F_r1-F1) + abs(Lm_r1-Lm1) + abs(Ds_r1-Ds) + abs(Ls_r1-Ls1) + abs(a_r1-a1) + abs(f_r1-f1) + abs(th_e_r1-th_e);
r2_ev = abs(Dm_r2-Dm) + abs(F_r2-F2) + abs(Lm_r2-Lm2) + abs(Ds_r2-Ds) + abs(Ls_r2-Ls2) + abs(a_r2-a2) + abs(f_r2-f2) + abs(th_e_r2-th_e);

if r1_ev < r2_ev
    a = a1;
    Lm = Lm1;
    f = f1;
    F = F1;
    Ls = Ls1;
else
    a = a2;
    Lm = Lm2;
    f = f2;
    F = F2;
    Ls = Ls2;
end
end