function a = ODReq22(Ls,e,beta,th_0)
a = (Ls)./( (2)+(((e.^2 - 1))/(e.*cos(beta - th_0) + 1)) );
end