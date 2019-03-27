function dFPR = ODReq10(h,Dm,f,beta)
dFPR = h - Dm./2 + 2.*f.*sin(beta);
end