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
% FF - FarField class 
% pol - {'x' | 'y' | 'lh' | 'rh'}
% th_M - subtended angle of the 'reflector' in rad

assert(strcmp(pol,'x')||strcmp(pol,'y')||strcmp(pol,'lh')||strcmp(pol,'rh')||isa('pol','double'),['Error: Unknown parameter for pol: ',pol])
% Get in BOR1
if ~FF.symmetryBOR1
    FF = FF.getBOR1pattern;
end
freq_vect = FF.freq;
lambda_vect = FF.c0./freq_vect;
k_vect = 2*pi./lambda_vect;

th = FF.th(1:FF.Ny);

[A1f,B1f,C1f,D1f] = FF.getBOR1comps;

[phi_0,phi_th_M,k_Delta,Delta,k_delta0,delta0,Z,eta_pd] = deal(zeros(1,FF.Nf));
for ff = 1:FF.Nf
    
    switch pol
        case 'x'
            CO = B1f(:,ff) + D1f(:,ff);
        case {'y','lh','rh'}
            CO = A1f(:,ff) + C1f(:,ff);
        otherwise
            error(['Unknown pol: ', pol])
    end
    
    % Move the pattern to the approximate PC (Kildal comments 1984)
    phi = unwrap(angle(CO));
    phi_0(ff) = phi(1);
    phi_th_M(ff) = interp1(th,phi,th_M);
    k_Delta(ff) = (phi_0(ff) - phi_th_M(ff))/(1 - cos(th_M));
    phi_Delta = phi - k_Delta(ff).*cos(th);
    Delta(ff) = k_Delta(ff)./k_vect(ff);
    
    % Find the PC from the formulas in Kildal 1983 (maximum eff method)
    % Weighting function
    w = abs(CO).*tan(th./2);
    % Integral constants (change number of points for very sharp patterns...)
    th_int = linspace(0,th_M,501);
    w_int = interp1(th,w,th_int);
    phi_Delta_int = interp1(th,phi_Delta,th_int);
    phi_Delta0 = interp1(th,phi_Delta,0);
    
    Iw = integral1D(th_int,w_int);
    Iwp = integral1D(th_int,w_int.*(phi_Delta_int - phi_Delta0));
    Iwp2 = integral1D(th_int,w_int.*(phi_Delta_int - phi_Delta0).^2);
    Iwc = integral1D(th_int,w_int.*(cos(th_int) - 1));
    Iwc2 = integral1D(th_int,w_int.*(cos(th_int) - 1).^2);
    Iwpc = integral1D(th_int,w_int.*(cos(th_int) - 1).*(phi_Delta_int - phi_Delta0));
    
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