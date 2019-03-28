function [f_root1,f_root2] = SDReq17(th_e,Ls,Ds)
%NEED TO FIND THE ROOT. THAT IS THE ANSWER.
%Equation: (4.*sin(th_e).*(4.*Ls.*tan(th_e./2) - Ds))Z^2
%+ (4.*Ds.*Ls(sin(th_e) + tan(th_e./2)) - 24.*Ls.^2.*tan(th_e./2).*sin(th_e))Z
%+ ((-Ls.^2.*Ds.*(2.*tan(th_e./2) + sin(th_e))) + Ls.^2.*tan(th_e./2).*sin(th_e).*(8.*Ls - tan(th_e./2).*Ds))

a = (4.*sin(th_e).*(4.*Ls.*tan(th_e./2) - Ds));
b = (4.*Ds.*Ls.*(sin(th_e) + tan(th_e./2)) - 24.*Ls.^2.*tan(th_e./2).*sin(th_e));
c = ((-Ls.^2.*Ds.*(2.*tan(th_e./2) + sin(th_e))) + Ls.^2.*tan(th_e./2).*sin(th_e).*(8.*Ls - tan(th_e./2).*Ds));

root1 = (-b + sqrt(b.^2 - 4.*a.*c))./(2.*a)
root2 = (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a)

f_root1 = root1;
f_root2 = root2;