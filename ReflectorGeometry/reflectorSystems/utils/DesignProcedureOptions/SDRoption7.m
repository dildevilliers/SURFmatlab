function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRoption7(Dm,Ds,Ls,th_e,sigma)
[F1,F2] = SDReq24(th_e,Ls,Ds,Dm,sigma);

f1 = SDReq15(th_e,Ls,F1,sigma,Dm);
a1 = SDReq14(th_e,Ls,F1,sigma,Dm);
Lm1 = SDReq12(Ls,F1,sigma,Dm,th_e);

f2 = SDReq15(th_e,Ls,F2,sigma,Dm);
a2 = SDReq14(th_e,Ls,F2,sigma,Dm);
Lm2 = SDReq12(Ls,F2,sigma,Dm,th_e);

[Dm_r1,F_r1,Lm_r1,Ds_r1,Ls_r1,a_r1,f_r1,th_e_r1] = SDRoption1(Dm,Lm1,Ls,th_e,sigma);
[Dm_r2,F_r2,Lm_r2,Ds_r2,Ls_r2,a_r2,f_r2,th_e_r2] = SDRoption1(Dm,Lm2,Ls,th_e,sigma);

r1_ev = abs(Dm_r1-Dm) + abs(F_r1-F1) + abs(Lm_r1-Lm1) + abs(Ds_r1-Ds) + abs(Ls_r1-Ls) + abs(a_r1-a1) + abs(f_r1-f1) + abs(th_e_r1-th_e);
r2_ev = abs(Dm_r2-Dm) + abs(F_r2-F2) + abs(Lm_r2-Lm2) + abs(Ds_r2-Ds) + abs(Ls_r2-Ls) + abs(a_r2-a2) + abs(f_r2-f2) + abs(th_e_r2-th_e);

if r1_ev < r2_ev
    a = a1;
    Lm = Lm1;
    f = f1;
    F = F1;
else
    a = a2;
    Lm = Lm2;
    f = f2;
    F = F2;
end

end