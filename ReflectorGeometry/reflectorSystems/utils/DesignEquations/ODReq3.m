function th_U = ODReq3(e,alpha,sigma,th_e,beta)
th_U = 2.*atan(((1 + e)./(1 - e)).*tan((alpha - sigma.*th_e)./2)) + beta;
end