function [F_root1,F_root2] = SDReq34(Df,Dm,sigma,Ls,th_e)
%NEED TO FIND THE ROOT. THAT IS THE ANSWER.
%Equation: (16.*tan(th_e).*(16.*Ls.^2.*(tan(th_e./2)).^2 + sigma.*Df.*Dm))Z^2
%+ (-8.*Dm.*(16.*Ls.^2.*sigma.*tan(th_e).*tan(th_e./2) + Df.*Dm))Z
%+ (tan(th_e).*Dm.^2.*(16.*Ls.^2 - sigma.*Df.*Dm))

a = (16.*tan(th_e).*(16.*Ls.^2.*(tan(th_e./2)).^2 + sigma.*Df.*Dm));
b = (-8.*Dm.*(16.*Ls.^2.*sigma.*tan(th_e).*tan(th_e./2) + Df.*Dm));
c = (tan(th_e).*Dm.^2.*(16.*Ls.^2 - sigma.*Df.*Dm));

root1 = (-b + sqrt(b.^2 - 4.*a.*c))./(2.*a);
root2 = (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a);

F_root1 = root1;
F_root2 = root2;