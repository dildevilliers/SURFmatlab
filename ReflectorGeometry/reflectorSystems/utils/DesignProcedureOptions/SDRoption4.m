function [Dm,F,Lm,Ds,Ls,a,f,th_e] = SDRoption4(F,Ds,Ls,th_e,sigma)
[Dm_root1,Dm_root2] = SDReq16(th_e,Ls,F,sigma,Ds);

Lm_1 = SDReq12(Ls,F,sigma,Dm_root1,th_e);
f_1 = SDReq15(th_e,Ls,F,sigma,Dm_root1);
a_1 = SDReq14(th_e,Ls,F,sigma,Dm_root1);

Lm_2 = SDReq12(Ls,F,sigma,Dm_root2,th_e);
f_2 = SDReq15(th_e,Ls,F,sigma,Dm_root2);
a_2 = SDReq14(th_e,Ls,F,sigma,Dm_root2);

[Dm_r1,F_r1,Lm_r1,Ds_r1,Ls_r1,a_r1,f_r1,th_e_r1] = SDRoption1(Dm_root1,Lm_1,Ls,th_e,sigma);
[Dm_r2,F_r2,Lm_r2,Ds_r2,Ls_r2,a_r2,f_r2,th_e_r2] = SDRoption1(Dm_root2,Lm_2,Ls,th_e,sigma);

r1_ev = abs(Dm_r1-Dm_root1) + abs(F_r1-F) + abs(Lm_r1-Lm_1) + abs(Ds_r1-Ds) + abs(Ls_r1-Ls) + abs(a_r1-a_1) + abs(f_r1-f_1) + abs(th_e_r1-th_e);
r2_ev = abs(Dm_r2-Dm_root2) + abs(F_r2-F) + abs(Lm_r2-Lm_2) + abs(Ds_r2-Ds) + abs(Ls_r2-Ls) + abs(a_r2-a_2) + abs(f_r2-f_2) + abs(th_e_r2-th_e);

if r1_ev < r2_ev
    Dm = Dm_root1;
    Lm = Lm_1;
    f = f_1;
    a = a_1;
else
    Dm = Dm_root2;
    Lm = Lm_2;
    f = f_2;
    a = a_2;
end
end