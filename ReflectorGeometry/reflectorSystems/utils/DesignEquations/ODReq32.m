function OP_0 = ODReq32(a,sigma,e,beta,th_0)
OP_0 = (-sigma.*a).*( (((e.^2 - 1))/(e.*cos(th_0 - beta) + 1)) );
end