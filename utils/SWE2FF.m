function FF = SWEgetField(Q,r,th,ph,F,P) %optional F matrix added!

% function FF = SWEgetField(Q,r,th,ph,P)
% Computes the farfield from a given set of spherical modes Q (as returned
% by readSPH.m) at distance r, and spherical angles th and ph.
% r, th, and ph must be the same lengths, and r can be inf for farfield
% requests
% P is an optional power (as a function of frequency) vector of size [1xNf]
% If it is provided it is simply put in the field structure FF as the
% radiated power.  P is also returned by readSPH.m

% Assume free space propagation
c0 = 299792458;
eta0 = 3.767303134749689e+02;    

if nargin > 5
    FF.Prad = P;
end

NF = length(Q.freq);

MMAX = [Q.MMAX];
NMAX = [Q.NMAX];
freq = [Q.freq];
if isempty(freq)
    warning('No frequency specified - working at 1 Hz');
    freq = ones(1,NF);
end
lam = c0./freq;
k = 2.*pi./lam;

% keyboard;

for ff = 1:NF
    disp(['Calculating farfield pattern ', num2str(ff),' of ', num2str(NF), '...'])
    
    Emat = F*Q.Q;
    
    
%     FF.Eth(:,ff) = eta0*(k(ff))*Emat(1,1,1,:,2);
%     FF.Eph(:,ff) = eta0*(k(ff))*Emat(1,1,1,:,3);
%     %BK: Removing scaling by free-space wave impedance
    FF.Eth(:,ff) = sqrt(eta0).*Emat(1:length(th),:);
    FF.Eph(:,ff) = sqrt(eta0).*Emat(length(th)+1:2*length(th),:);
end

FF.th = repmat(th,1,NF);
FF.ph = repmat(ph,1,NF);
if ~all(isinf(r)) 
    FF.r = repmat(r,1,NF);
end
if ~isempty(freq)
    FF.freq = freq;
end
