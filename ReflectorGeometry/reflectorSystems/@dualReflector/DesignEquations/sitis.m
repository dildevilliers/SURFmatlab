function f = sitis(dFPR,h,Dm,beta)
f = (dFPR - h + DM./2)./(2.*sin(beta));
end