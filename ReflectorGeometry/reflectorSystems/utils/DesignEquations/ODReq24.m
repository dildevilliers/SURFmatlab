function h = ODReq24(dFPR,Dm,f,beta)
h = dFPR + Dm./2 - 2.*f.*sin(beta);
end