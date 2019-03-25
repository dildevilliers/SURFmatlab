function Ds = SDReq7(Ls,f,th_e,sigma,F,Dm)
Ds = (4.*(Ls - f))./( (1)./(sin(th_e)) + (sigma.*(16.*F.^2 + Dm.^2))./(8.*F.*Dm) )
end