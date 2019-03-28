function OP_1 = ODReq33(a,sigma,e,beta,th_L)
OP_1 = (-sigma.*a).*( (((e.^2 - 1))/(e.*cos(th_L - beta) + 1)) );
end