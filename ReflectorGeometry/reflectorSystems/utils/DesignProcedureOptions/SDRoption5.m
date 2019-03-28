function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRoption5(Lm,Ds,Ls,th_e,sigma)
[f1,f2] = SDReq17(th_e,Ls,Ds);

Dm1 = SDReq18(th_e,Lm,sigma,Ls,f1);
F1 = SDReq5(Lm,f1);
a1 = SDReq6(Ls,f1);

Dm2 = SDReq18(th_e,Lm,sigma,Ls,f2);
F2 = SDReq5(Lm,f2);
a2 = SDReq6(Ls,f2);

[Dm_r1,F_r1,Lm_r1,Ds_r1,Ls_r1,a_r1,f_r1,th_e_r1] = SDRoption1(Dm1,Lm,Ls,th_e,sigma);
[Dm_r2,F_r2,Lm_r2,Ds_r2,Ls_r2,a_r2,f_r2,th_e_r2] = SDRoption1(Dm2,Lm,Ls,th_e,sigma);

r1_ev = abs(Dm_r1-Dm1) + abs(F_r1-F1) + abs(Lm_r1-Lm) + abs(Ds_r1-Ds) + abs(Ls_r1-Ls) + abs(a_r1-a1) + abs(f_r1-f1) + abs(th_e_r1-th_e);
r2_ev = abs(Dm_r2-Dm2) + abs(F_r2-F2) + abs(Lm_r2-Lm) + abs(Ds_r2-Ds) + abs(Ls_r2-Ls) + abs(a_r2-a2) + abs(f_r2-f2) + abs(th_e_r2-th_e);

if r1_ev < r2_ev
    Dm = Dm1;
    F = F1;
    f = f1;
    a = a1;
else
    Dm = Dm2;
    F = F2;
    f = f2;
    a = a2;
end

end