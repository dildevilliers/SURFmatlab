function [F_root1,F_root2] = SDReq24(th_e,Ls,Ds,Dm,sigma)
%NEED TO FIND THE ROOT. THAT IS THE ANSWER.
%Equation: (16.*sin(th_e).*(4.*Ls.*tan(th_e./2) - Ds))Z^2
%+ (8.*Dm.*sigma.*(2.*sin(th_e).*Ls - Ds))Z
%+ (-Ds.*Dm.^2.*sin(th_e))

a = (16.*sin(th_e).*(4.*Ls.*tan(th_e./2) - Ds));
b = (8.*Dm.*sigma.*(2.*sin(th_e).*Ls - Ds));
c = (-Ds.*Dm.^2.*sin(th_e));

root1 = (-b + sqrt(b.^2 - 4.*a.*c))./(2.*a);
root2 = (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a);

F_root1 = root1;
F_root2 = root2;