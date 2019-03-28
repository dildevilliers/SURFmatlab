function f = SDReq15(th_e,Ls,F,sigma,Dm)
f = (Ls.*(sigma.*Dm - 4.*F.*tan(th_e./2)))./(2.*sigma.*Dm);
end