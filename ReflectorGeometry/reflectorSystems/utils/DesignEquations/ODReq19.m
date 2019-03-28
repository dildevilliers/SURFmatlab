function a = ODReq19(Lt,h,Dm,F,sigma,e,beta,th_U,th_L)
a = (Lt + (((2.*h - Dm).^2)/(16.*F)) - F)./( (((sigma - 1)/2).*((e.^2 - 1).*(cos(th_U)))/(e.*cos(beta - th_U) + 1))-(((sigma + 1)/2).*((e.^2 - 1).*(cos(th_L)))/(e.*cos(beta - th_L) + 1)) );
end