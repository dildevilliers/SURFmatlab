function a = ODReq20(Ht,h,Dm,sigma,e,beta,th_U,th_L)
a = (Ht - h - Dm./2)./( (((sigma + 1)/2).*((e.^2 - 1).*(sin(th_U)))/(e.*cos(beta - th_U) + 1))-(((sigma - 1)/2).*((e.^2 - 1).*(sin(th_L)))/(e.*cos(beta - th_L) + 1)) );
end