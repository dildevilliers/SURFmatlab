function Ds = SDReq13(th_e,Ls,F,sigma,Dm)
Ds = (16.*sin(th_e).*F.*Ls.*(sigma.*Dm + 4.*F.*tan(th_e./2)))./(8.*sigma.*F.*Dm + sin(th_e).*(Dm.^2 + 16.*F.^2));
end