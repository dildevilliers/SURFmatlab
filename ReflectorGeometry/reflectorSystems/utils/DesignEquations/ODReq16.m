function f = ODReq16(dFPR,h,Dm,beta)
f = (dFPR - h + Dm./2)./(2.*sin(beta));
end