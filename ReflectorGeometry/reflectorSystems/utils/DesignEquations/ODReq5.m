function e = ODReq5(sigma,beta,th_0)
e = (1 - sigma.*sqrt((tan(beta./2))./(tan((beta - th_0)./2))))./(1 + sigma.*sqrt((tan(beta./2))./(tan((beta - th_0)./2))));
end