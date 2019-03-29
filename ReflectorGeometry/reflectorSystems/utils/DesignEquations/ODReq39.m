function Dsy = ODReq39(a,e,th_e,alpha)
ph = linspace(0,2*pi,1001);
Dsy = max( (2.*a.*(e.^2 - 1).*sin(th_e).*sin(ph))./(e.*(-sin(alpha).*sin(th_e).*cos(ph) + cos(alpha).*cos(th_e)) - 1) );
end