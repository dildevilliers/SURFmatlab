function th_e = SDReq26(F,Dm,Ds,f,sigma)
th_e = atan( (8.*F.*Dm.*Ds)./(32.*f.*F.*Dm + sigma.*Ds.*(16.*F.^2 - Dm.^2)) );
end