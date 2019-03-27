function [Dm_root1,Dm_root2] = SDReq16(th_e,Ls,F,sigma,Ds)
%NEED TO FIND THE ROOT. THAT IS THE ANSWER.
%Equation: (Ds.*sin(th_e))Z^2 + (8.*F.*sigma(Ds - 2.*sin(th_e).*Ls))Z + (16.*F.^2.*sin(th_e).*(Ds - 4.*Ls.*tan(th_e./2)))
a = (Ds.*sin(th_e));
b = (8.*F.*sigma.*(Ds - 2.*sin(th_e).*Ls));
c = (16.*F.^2.*sin(th_e).*(Ds - 4.*Ls.*tan(th_e./2)));
root1 = (-b + sqrt(b.^2 - 4.*a.*c))./(2.*a);
root2 = (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a);

Dm_root1 = root1;
Dm_root2 = root2;
