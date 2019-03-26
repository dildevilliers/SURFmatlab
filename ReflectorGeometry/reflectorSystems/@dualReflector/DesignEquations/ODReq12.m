function Lt = ODReq12(a,sigma,e,th_L,beta,th_U,Dm,F,h)
Lt = -a.*((sigma + 1)./2).*(((e.^2 - 1).*cos(th_L))./(e.*cos(beta - th_L) + 1)) + a.*((sigma - 1)./2).*(((e.^2 - 1).*cos(th_U))./(e.*cos(beta - th_U) + 1)) - ((2.*h - Dm).^2)./(16.*F) + F;
end