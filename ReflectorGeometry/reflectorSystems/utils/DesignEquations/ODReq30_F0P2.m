function F0P2 = ODReq30_F0P2(a,sigma,e,beta,th_U)
OP_2 = ODReq34(a,sigma,e,beta,th_U);
F0P2 = 2.*a - sigma.*OP_2;
end