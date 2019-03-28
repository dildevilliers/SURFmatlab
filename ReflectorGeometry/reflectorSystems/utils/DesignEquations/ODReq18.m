function a = ODReq18(sigma,Dsx,e,beta,th_U,th_L)
a = (-sigma.*Dsx)./( (((e.^2 - 1).*(sin(beta - th_U)))/(e.*cos(beta - th_U) + 1))-(((e.^2 - 1).*(sin(beta - th_L)))/(e.*cos(beta - th_L) + 1)) );
end