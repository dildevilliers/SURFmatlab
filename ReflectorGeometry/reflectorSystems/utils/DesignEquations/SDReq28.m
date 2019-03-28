function f = SDReq28(F,Dm,Df,sigma,th_e)
f = sqrt( (8.*F.*Dm.*Df - sigma.*Df.*tan(th_e).*(16.*F.^2 - Dm.^2))./(64.*tan(th_e).*Dm) );
end