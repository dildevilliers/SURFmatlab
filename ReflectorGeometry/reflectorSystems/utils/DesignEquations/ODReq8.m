function Ls = ODReq8(a,e,beta,th_0)
Ls = a.*(2 + ((e.^2 - 1)./(e.*cos(beta - th_0) + 1)));
end