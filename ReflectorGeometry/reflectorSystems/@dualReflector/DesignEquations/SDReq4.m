function f = SDReq4(Ls,sigma,Dm,Lm,th_e)
f = Ls.*( (sigma.*Dm - 4.*Lm.*tan(th_e./2))./(2.*sigma.*Dm + 8.*Ls.*tan(th_e./2)) );
end