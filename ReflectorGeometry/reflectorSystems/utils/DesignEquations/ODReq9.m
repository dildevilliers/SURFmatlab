function Lm = ODReq9(a,e,beta,th_0,h)
Lm = -a.*(((e.^2 - 1)./(e.*cos(beta - th_0) + 1))) - h./sin(th_0);
end