function [F_root1,F_root2] = SDReq32(Df,Dm,sigma,Ds,th_e)
%NEED TO FIND THE ROOT. THAT IS THE ANSWER.
%Equation: (16.*(Df.*Dm + sigma.*Ds.^2))Z^2
%+ ((-8.*Dm.*Ds.^2)./(tan(th_e)))Z
%+ (-sigma.*Ds.^2.*Dm.^2)

a = (16.*(Df.*Dm + sigma.*Ds.^2));
b = ((-8.*Dm.*Ds.^2)./(tan(th_e)));
c = (-sigma.*Ds.^2.*Dm.^2);

root1 = (-b + sqrt(b.^2 - 4.*a.*c))./(2.*a);
root2 = (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a);

F_root1 = root1;
F_root2 = root2;