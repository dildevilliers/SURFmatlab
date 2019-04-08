function Dsx = ODReq27(a,sigma,e,beta,th_U,th_L)
Dsx = (-sigma.*a).*( (((e.^2 - 1).*(sin(beta - th_U)))/(e.*cos(beta - th_U) + 1))-(((e.^2 - 1).*(sin(beta - th_L)))/(e.*cos(beta - th_L) + 1)) );
end