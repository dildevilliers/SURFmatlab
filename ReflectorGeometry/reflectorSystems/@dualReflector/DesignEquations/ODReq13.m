function Ht = ODReq13(h,Dm,a,sigma,e,th_L,beta,th_U)
Ht = h + Dm./2 - a.*((sigma - 1)./2).*(((e.^2 - 1).*sin(th_L))./(e.*cos(beta - th_L) + 1)) + a.*((sigma + 1)./2).*(((e.^2 - 1).*sin(th_U))./(e.*cos(beta - th_U) + 1));
end