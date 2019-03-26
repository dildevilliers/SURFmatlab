function f = SDReq35(Ls,th_e,Dm,F,sigma)
%IMPORTANT NOTE: THE PAPER SPECIFIES TH_2 BUT I THINK IT SHOULD BE TH_E.
%TH_E IS THUS IMPLEMENTED EVERYWHERE IN THIS EQUATION INSTEAD OF TH_2.
f = ((-1/2).*Ls).*( (4.*F.*tan(th_e./2) - sigma.*Dm)./(sigma.*Dm) );
end