function dP = dPcos_dthN(th,m,n)

% Calculates the analytical derivative of the Normalised Legendre function
% used for SWE: dP = d/dth(P^m_n(cos(th)))

% Simplified form - from chain rule
% Notice the normalisation in the second term below.  Do like this to avoid
% large number factorial divisions and trust internal matlab handles this
% well for legendre.m (not sure it does...)

dP = -csc(th).*((n+1).*cos(th).*PmnNorm(cos(th),m,n) + (m-n-1).*PmnNorm(cos(th),m,n+1)/sqrt((n+1+0.5)/(n+0.5)*(n+1-m)/(n+1+m)));

% Sort out indeterminate forms
iNaN0 = find(th == 0); 
iNaNpi = find(th == pi); 
if m == 1
    dP(iNaN0) = n.*(n+1)/2*sqrt((n+0.5)/(n*(n+1))); % McLauren expansion (wolframalpha)
    dP(iNaNpi) = (-1)^n.*dP(iNaN0(1));
else
    dP([iNaN0,iNaNpi]) = 0;
end
