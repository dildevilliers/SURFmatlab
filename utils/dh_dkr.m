function dh = dh_dkr(x,n)

% Calculates the analytical derivative of the Spherical Hankel function
% used for SWE: dh = 1/x*(d/dx(xh^(2)_n(x)))

if n == 0
    dh = exp(-1i.*x);
else
    % Recurrence relation
    dh = h2Sph(x,n-1) - n.*h2Sph(x,n)./x;
end
    
%% Farfield asymptote
dh(isinf(x)) = 1i^n; % Suppress the e^(-jkr)/kr factor


