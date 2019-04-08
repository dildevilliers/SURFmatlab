function Dm = SDReq18(th_e,Lm,sigma,Ls,f)
Dm = (4.*Ls.*tan(th_e./2).*(Lm + 2.*f))./(sigma.*(Ls - 2.*f));
end