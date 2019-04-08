function dP = dPcos_dth(th,m,n)

% Calculates the analytical derivative of the Un-normalised Legendre function
% used for SWE: dP = d/dth(P^m_n(cos(th)))

% Straight chain rule
% dP = -1./(cos(th).^2 - 1).*(sin(th).*((-n-1).*cos(th).*Pmn(cos(th),m,n) + (-m+n+1).*Pmn(cos(th),m,n+1)));
% Simplified form
dP = -csc(th).*((n+1).*cos(th).*Pmn(cos(th),m,n) + (m-n-1).*Pmn(cos(th),m,n+1));
% Sort out indeterminate forms
iNaN0 = find(th == 0); 
iNaNpi = find(th == pi); 
if m == 1
    dP(iNaN0) = -n.*(n+1)/2; % McLauren expansion (wolframalpha)
    dP(iNaNpi) = (-1)^n.*dP(iNaN0(1));
else
    dP([iNaN0,iNaNpi]) = 0;
end

