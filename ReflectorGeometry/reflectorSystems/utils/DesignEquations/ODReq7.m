function th_e = ODReq7(sigma,e,th_U,beta,alpha)
th_e = -sigma.*(2.*atan(((1 - e)./(1 + e)).*tan((th_U - beta)./2)) - alpha);
end