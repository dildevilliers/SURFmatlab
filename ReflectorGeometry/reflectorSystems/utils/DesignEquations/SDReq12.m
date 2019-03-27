function Lm = SDReq12(Ls,F,sigma,Dm,th_e)
Lm =  -((sigma.*Dm).*(Ls-F) - 4.*Ls.*F.*tan(th_e./2))./(sigma.*Dm);
end