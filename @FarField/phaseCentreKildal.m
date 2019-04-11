function [Z, Delta, delta0, eta_pd] = phaseCentreKildal(FF,pol,th_M)

% function [PC, eta_phi] = PCoptFunc(pol,freq_vect,th_M,BOR1pattern)
% computes the phase center and approximate phase efficiency of a given
% farfield pattern using Kildal 1983 paper
% Returns:
% Z - the total phase center (Delta + delta0) position in [m] as a frequency vector
% Delta - the initial approximation of the PC [m]
% delta0 - the fine approximation around Delta of the PC [m]
% eta_pd - the phase efficiency pu as a frequency vector (NaN for large PC displacements)
% Inputs:
% pol - 'x' or 'y'
% freq_vect - in Hz
% th_M - subtended angle of the 'reflector' in rad
% BOR1pattern - structure containing the variables th (vector) and A1 and C1 - each a
% matrix of size [nr_th, nr_freq].  0 <= th <= pi.

load constants

th = BOR1pattern.th;

% Rough error checking...
if (th(1) ~= 0) | (th(length(th)) ~= pi)
    error('Theta should range from 0 to pi!');
end

% Set up propagation constant
lambda_vect = c0./freq_vect;
k_vect = 2*pi./lambda_vect;
Nf = length(freq_vect);

[phi_0,phi_th_M,k_Delta,Delta,k_delta0,delta0,Z,eta_pd] = deal(zeros(1,Nf));
for ff = 1:Nf
    
    if pol == 'x'
        A1 = BOR1pattern.B1(:,ff);
        C1 = BOR1pattern.D1(:,ff);
    elseif pol == 'y'
        A1 = BOR1pattern.A1(:,ff);
        C1 = BOR1pattern.C1(:,ff);
    end
    
    % Set up th vector for the current frequency
    th_calc = th(th <= th_M);
    
    % Move the pattern to the approximate PC (Kildal comments 1984)
    phi = unwrap(angle(A1 + C1));
    phi_0(ff) = phi(1);
    phi_th_M(ff) = interp1(th,phi,th_M);
    k_Delta(ff) = (phi_0(ff) - phi_th_M(ff))/(1 - cos(th_M));
    phi_Delta = phi - k_Delta(ff).*cos(th);
    Delta(ff) = k_Delta(ff)./k_vect(ff);
    
    % Find the PC from the formulas in Kildal 1983 (maximum eff method)
    % Weighting function
    w = abs(A1 + C1).*tan(th./2);
    % Integral constants (change number of points for very sharp patterns...)
    th_int = linspace(0,th_M,501);
    dth = th_int(2) - th_int(1);
    w_int = interp1(th,w,th_int);
    phi_Delta_int = interp1(th,phi_Delta,th_int);
    phi_Delta0 = interp1(th,phi_Delta,0);
    
    Iw = trapz(w_int).*dth;
    Iwp = trapz(w_int.*(phi_Delta_int - phi_Delta0)).*dth;
    Iwp2 = trapz(w_int.*(phi_Delta_int - phi_Delta0).^2).*dth;
    Iwc = trapz(w_int.*(cos(th_int) - 1)).*dth;
    Iwc2 = trapz(w_int.*(cos(th_int) - 1).^2).*dth;
    Iwpc = trapz(w_int.*(cos(th_int) - 1).*(phi_Delta_int - phi_Delta0)).*dth;
    
    % PC pos
    k_delta0(ff) = (Iw.*Iwpc - Iwp.*Iwc)./(Iwc2.*Iw - Iwc.^2);
    delta0(ff) = k_delta0(ff)./k_vect(ff);
    
    Z(ff) = Delta(ff) + delta0(ff);
    
    % Calculate the approximate phase efficiency
    if abs(Z(ff)/lambda_vect(ff)) > pi/4   % Return NaN for large PC errors - approximation not valid
        eta_pd(ff) = NaN;
    else
        c = 1 - Iwp2/Iw + (Iwp/Iw)^2;
        b = Iwpc/Iw - (Iwp*Iwc/Iw^2);
        a = Iwc2/Iw - (Iwc/Iw)^2;
        
        eta_pd(ff) = c + 2*b*k_vect(ff)*Z(ff) - a*(k_vect(ff)*Z(ff)).^2;
    end
end