function F0P1 = ODReq30_F0P1(a,sigma,e,beta,th_L)
OP_1 = ODReq33(a,sigma,e,beta,th_L);
F0P1 = 2.*a - sigma.*OP_1;
end