function [F_root1,F_root2] = SDReq33(Df,Dm,sigma,Lm,th_e)
%NEED TO FIND THE ROOT. THAT IS THE ANSWER.
%Equation: (16.*tan(th_e).*(Dm + sigma.*Df))Z^2
%+ (-8.*Dm.*(4.*Lm.*tan(th_e) + Df))Z
%+ (tan(th_e).*Dm.*(16.*Lm.^2 - sigma.*Df.*Dm))

a = (16.*tan(th_e).*(Dm + sigma.*Df));
b = (-8.*Dm.*(4.*Lm.*tan(th_e) + Df));
c = (tan(th_e).*Dm.*(16.*Lm.^2 - sigma.*Df.*Dm));

root1 = (-b + sqrt(b.^2 - 4.*a.*c))./(2.*a);
root2 = (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a);

F_root1 = root1;
F_root2 = root2;