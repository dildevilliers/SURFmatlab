function F = SDReq1(Ls,sigma,Dm,Lm,th_e)
F = Ls.*( (sigma.*Dm - 4.*Lm.*tan(th_e./2))./(2.*sigma.*Dm + 8.*Ls.*tan(th_e./2)) );
end