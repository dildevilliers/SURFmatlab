function Ls = SDReq9(Lm,F,sigma,Dm,th_e)
Ls =  -((sigma.*Dm).*(Lm-F))./(sigma.*Dm - 4.*F.*tan(th_e./2));
end