function Dm = SDReq16(th_e,Ls,F,sigma,Ds)
%NEED TO FIND THE ROOT. THAT IS THE ANSWER.
%Equation: (Ds.*sin(th_e))Z^2 + (8.*F.*sigma(Ds - 2.*sin(th_e).*Ls))Z + (16.*F.^2.*sin(th_e).*(Ds - 4.*Ls.*tan(th_e./2)))
a = (Ds.*sin(th_e));
b = (8.*F.*sigma(Ds - 2.*sin(th_e).*Ls));
c = (16.*F.^2.*sin(th_e).*(Ds - 4.*Ls.*tan(th_e./2)));
root1 = (-b + sqrt(b.^2 - 4.*a.*c))./(2.*a);
root2 = (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a);

%This next part is improvised. I'm not sure which root to take.
if root1 >= 0
    Dm = root1;
else
    Dm = root2;
end