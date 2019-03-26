function OP_2 = ODReq34(a,sigma,e,beta,th_U)
OP_2 = (-sigma.*a).*( (((e.^2 - 1))/(e.*cos(th_U - beta) + 1)) );
end