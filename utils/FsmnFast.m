function [Fsm0n,Fsmmn,Fsmpn] = FsmnFast(r,th,ph,MMAX,NMAX,f)

% [Fsm0n,Fsmmn,Fsmpn] = FsmnFast(r,th,ph,MMAX,NMAX,f)
% Function to compute the F factors for the Q-SWE in (3.129) and (3.130) in
% the GRASP Technical Description
% Inputs:
% r,th,ph: Spherical coordinates [m,rad,rad].  Can be vectors of the same
%          lengths, and r = inf implies a farfield limit with the
%          exponential phase/amplitude suppressed [e^{jkr}/kr]
% MMAX,NMAX: Maximum mode number integers
% f: Scalar frequency in [Hz]
% Outputs:
% Fsm0n,Fsmmn,Fsmpn: [2,MMAX,NMAX,Ndirections,3] matrix with last 'columns' the [r,th,ph] components of the vector
%                       function.  All the modes are output for each call.

% Basic error checking
size_r = size(r);
size_th = size(th);
size_ph = size(ph);
assert(isequal(size_r,size_th) & isequal(size_th,size_ph),'Size of r, th and ph must be the same');
assert(length(f) == 1,'Frequency (f) should be a scalar')
NDIR = prod(size_r);

% Force column vectors
r = reshape(r,1,1,NDIR);
th = reshape(th,1,1,NDIR);
ph = reshape(ph,1,1,NDIR);
% r = r(:);
% th = th(:);
% ph = ph(:);

% Build the MMAX x NMAX x NDIR arrays
r = repmat(r,MMAX,NMAX,1);
th = repmat(th,MMAX,NMAX,1);
ph = repmat(ph,MMAX,NMAX,1);

[Nmat,Mmat] = meshgrid(1:NMAX,1:MMAX);
Nmat = repmat(Nmat,1,1,NDIR);
Mmat = repmat(Mmat,1,1,NDIR);

% Initial admin
l = 299792458/f;
k = 2*pi/l;
kr = k.*r;

% Get m_mn and n_mn quickly with no recalculation of special functions
% Use recursion as far as possible and put functions inline which uses
% other versions of the same special function

% Initialize
[hkr,dhkr] = deal(zeros(MMAX,NMAX,NDIR));
[Pcell] = deal(cell(1,NMAX));
% Get zero/first order cases for recurrence
ekr = exp(-1i.*kr)./kr;
ekr(isinf(kr)) = 1; % Suppress for farfield calculations
% h0kr(:,1,:) = 1i.*ekr(:,1,:);
h0kr = 1i.*ekr(:,1,:);
hkr(:,1,:) = (-1 + 1i./kr(:,1,:)).*ekr(:,1,:);
% keyboard;
dhkr(:,1,:) = h0kr(:,1,:) - hkr(:,1,:)./kr(:,1,:);
costh_vect = reshape(cos(th(1,1,:)),NDIR,1);
Pcell{1} = legendre(1,costh_vect,'norm');
for nn = 2:NMAX
    % Hankel
    % Recurrence relation
    if nn == 2
        hkr(:,nn,:) = (2*(nn-1) + 1).*hkr(:,nn-1,:)./kr(:,nn,:) - h0kr(:,1,:);
    else
        hkr(:,nn,:) = (2*(nn-1) + 1).*hkr(:,nn-1,:)./kr(:,nn,:) - hkr(:,nn-2,:);
    end
    % Farfield asymptote
    hkr(:,nn,isinf(kr(1,nn,:))) = 1i^(nn+1); % Suppress the e^(-jkr)/kr factor
    % Hankel derivative
    dhkr(:,nn,:) = hkr(:,nn-1,:) - nn.*hkr(:,nn,:)./kr(:,nn,:);
    % Farfield asymptote
    dhkr(:,nn,isinf(kr(1,nn,:))) = 1i^nn; % Suppress the e^(-jkr)/kr factor
    
    % Legendre - Get for all m <= n since it is provided for free anyway...
    Pcell{nn} = legendre(nn,costh_vect,'norm');
end

% Repack the Legendre cells and get derivatives
[Pcosth,dPcosth] = deal(zeros(MMAX,NMAX,NDIR));
[Pcosth0,dPcosth0] = deal(zeros(1,NMAX,NDIR));
% Sort out indeterminate forms
iNaN0 = find(th(1,1,:) == 0);
iNaNpi = find(th(1,1,:) == pi);
cscth_vect = reshape(csc(th(1,1,:)),NDIR,1);

Min1Fact = (-1); % Should be -1 for GRASP and 1 for FEKO

for nn = 1:NMAX
    for mm = 1:min(nn,MMAX)
        
        if nn == 1  % Only calculate on n = 1
            Pcosth0(1,1,:) = Pcell{1}(1,:);
            Pcosth(mm,1,:) = Pcell{1}(mm+1,:);
        end
        
        Pcosth0(1,nn,:) = Pcell{nn}(1,:);
        Pcosth(mm,nn,:) = Pcell{nn}(mm+1,:);
        if nn < NMAX
            Pcosth0(1,nn+1,:) = Pcell{nn+1}(1,:);
            Pcosth(mm,nn+1,:) = Pcell{nn+1}(mm+1,:);
            
%             -(Min1Fact).^mm.*
            
            if mm <= nn
%                 dPcosth0(1,nn,:) = -cscth_vect.*((nn+1).*costh_vect.*reshape(Pcosth0(1,nn,:),NDIR,1) + (nn-1).*reshape(Pcosth0(1,nn+1,:),NDIR,1)/sqrt((nn+1+0.5)/(nn+0.5)*(nn+1)/(nn+1)));
                dPcosth0(1,nn,:) = -sqrt(nn+0.5).*(nn+1).*cscth_vect.*(costh_vect.*reshape(Pcosth0(1,nn,:),NDIR,1)./sqrt(nn+0.5) - reshape(Pcosth0(1,nn+1,:),NDIR,1)./sqrt(nn+1.5));
                dPcosth(mm,nn,:) = -cscth_vect.*((nn+1).*costh_vect.*reshape(Pcosth(mm,nn,:),NDIR,1) + (mm-nn-1).*reshape(Pcosth(mm,nn+1,:),NDIR,1)/sqrt((nn+1+0.5)/(nn+0.5)*(nn+1-mm)/(nn+1+mm)));

                % Handle undetermined cases
                if mm == 1
                    if ~isempty(iNaN0)
                        dPcosth(mm,nn,iNaN0) = (nn).*(nn+1)/2*sqrt((nn+0.5)/((nn)*(nn+1))); % McLauren expansion (wolframalpha)
                        if ~isempty(iNaNpi)
                            dPcosth(mm,nn,iNaNpi) = (-1)^(nn).*dPcosth(mm,nn,iNaN0(1));
                        end
                    end
                else
                    dPcosth(mm,nn,[iNaN0,iNaNpi]) = 0;
                end
                dPcosth0(1,nn,[iNaN0,iNaNpi]) = 0;
            end
            % Get the last case where the recursion lookup table runs out
        else
%             Pcosth0(1,nn,:) = Pcell{nn}(1,:);
%             Pcosth(mm,nn,:) = (Min1Fact).^mm.*Pcell{nn}(mm+1,:);
%             
            PcosthNp1mat = legendre(nn+1,costh_vect,'norm');
            PcosthNp1 = reshape(PcosthNp1mat(mm+1,:),NDIR,1);
            Pcosth0Np1 = reshape(PcosthNp1mat(1,:),NDIR,1);
            %                 keyboard;
%             dPcosth0(1,nn,:) = -cscth_vect.*((nn+1).*costh_vect.*reshape(Pcosth0(1,nn,:),NDIR,1) + (nn-1).*Pcosth0Np1/sqrt((nn+1.5)/(nn+0.5)*(nn+1)/(nn+1)));
            dPcosth0(1,nn,:) = -sqrt(nn+0.5).*(nn+1).*cscth_vect.*(costh_vect.*reshape(Pcosth0(1,nn,:),NDIR,1)./sqrt(nn+0.5) - Pcosth0Np1./sqrt(nn+1.5));
            dPcosth(mm,nn,:) = -cscth_vect.*((nn+1).*costh_vect.*reshape(Pcosth(mm,nn,:),NDIR,1) + (mm-nn-1).*PcosthNp1/sqrt((nn+1.5)/(nn+0.5)*(nn+1-mm)/(nn+1+mm)));
            % Handle undetermined cases
            if mm == 1
                dPcosth(mm,nn,iNaN0) = nn.*(nn+1)/2*sqrt((nn+0.5)/(nn*(nn+1))); % McLauren expansion (wolframalpha)
%                 dPcosth(mm,nn,iNaNpi) = (-1)^nn.*dPcosth(mm,nn,iNaN0(1));
                dPcosth(mm,nn,iNaNpi) = (-1)^nn.*nn.*(nn+1)/2*sqrt((nn+0.5)/(nn*(nn+1)));
            else
                dPcosth(mm,nn,[iNaN0,iNaNpi]) = 0;
            end
            dPcosth0(1,nn,[iNaN0,iNaNpi]) = 0;
        end
        
%         % Test using functions...
%         P0test = PmnNorm(reshape(cos(th(1,nn,:)),NDIR,1),0,nn);
%         P0err(mm,nn) = max(abs(P0test - reshape(Pcosth0(1,nn,:),NDIR,1)))
%         
%         dP0test = dPcos_dthN(reshape(th(1,nn,:),NDIR,1),0,nn);
%         [dP0err(mm,nn), idP0errMax(mm,nn)] = max(abs(dP0test - reshape(dPcosth0(1,nn,:),NDIR,1)))
%         
%         if mm <= nn
%             Ptest = PmnNorm(reshape(cos(th(mm,nn,:)),NDIR,1),mm,nn);
%             Perr(mm,nn) = max(abs(Ptest - reshape(Pcosth(mm,nn,:),NDIR,1)))
%             figure(1)
%             plot(Ptest,'b'), grid on, hold on
%             plot(reshape(Pcosth(mm,nn,:),NDIR,1),'k-')
%             hold off
%             
%             dPtest = dPcos_dthN(reshape(th(mm,nn,:),NDIR,1),mm,nn);
%             [dPerr(mm,nn), idPerrMax(mm,nn)] = max(abs(dPtest - reshape(dPcosth(mm,nn,:),NDIR,1)))
%             figure(2)
%             plot(dPtest,'b'), grid on, hold on
%             plot(reshape(dPcosth(mm,nn,:),NDIR,1),'k-')
%             hold off
%             
%         end
%         keyboard;
%         
    end
end

% keyboard

% Save some memory
clear Pcell

%% Get the mprime and nprimes (eq 3.131 and 3.132 in GRASP technical description)

ephm = exp(Min1Fact.*1i.*(-Mmat).*ph);
ephp = exp(Min1Fact.*1i.*(+Mmat).*ph);


% mprime
% r-component: All zeros so handled later to save memory
% th-component
% GRASP Manual signs
% mp_mmn_th_num = -hkr.*1i.*(-Mmat).*Pcosth.*ephm;    %-m
% mp_mpn_th_num = -hkr.*1i.*(+Mmat).*Pcosth.*ephp;    %+m
% Hansen/FEKO signs
mp_mmn_th_num = Min1Fact.*hkr.*1i.*(-Mmat).*Pcosth.*ephm;    %-m
mp_mpn_th_num = Min1Fact.*hkr.*1i.*(+Mmat).*Pcosth.*ephp;    %+m

% No need for m=0 case - all zeros and handled later
mp_mmn_th = mp_mmn_th_num./sin(th);
mp_mpn_th = mp_mpn_th_num./sin(th);
% Sort out indeterminate forms by L'Hospital
iNaNm = find(abs(mp_mmn_th_num) < eps & abs(sin(th)) < eps);
iNaNp = find(abs(mp_mpn_th_num) < eps & abs(sin(th)) < eps);
% GRASP manual signs
% mp_mmn_th(iNaNm) = -hkr(iNaNm).*1i.*(-Mmat(iNaNm)).*dPcosth(iNaNm).*ephm(iNaNm)./cos(th(iNaNm));
% mp_mpn_th(iNaNp) = -hkr(iNaNp).*1i.*(+Mmat(iNaNp)).*dPcosth(iNaNp).*ephp(iNaNp)./cos(th(iNaNp));
% Hansen/FEKO signs
mp_mmn_th(iNaNm) = Min1Fact.*hkr(iNaNm).*1i.*(-Mmat(iNaNm)).*dPcosth(iNaNm).*ephm(iNaNm)./cos(th(iNaNm));
mp_mpn_th(iNaNp) = Min1Fact.*hkr(iNaNp).*1i.*(+Mmat(iNaNp)).*dPcosth(iNaNp).*ephp(iNaNp)./cos(th(iNaNp));
% keyboard;
clear mp_mmn_th_num mp_mpn_th_num
clear iNaNm iNaNp
% ph-component
mp_mmn_ph = -hkr.*dPcosth.*ephm;    %-m
mp_mpn_ph = -hkr.*dPcosth.*ephp;    %+m
mp_m0n_ph = -hkr(1,:,:).*dPcosth0;    %m=0


% nprime
% r-component
np_mmn_r = Nmat.*(Nmat+1)./kr.*hkr.*Pcosth.*ephm;    %-m
np_mpn_r = Nmat.*(Nmat+1)./kr.*hkr.*Pcosth.*ephp;    %+m
np_m0n_r = Nmat(1,:,:).*(Nmat(1,:,:)+1)./kr(1,:,:).*hkr(1,:,:).*Pcosth0;    %m=0

% th-component
np_mmn_th = dhkr.*dPcosth.*ephm;    %-m
np_mpn_th = dhkr.*dPcosth.*ephp;    %+m
np_m0n_th = dhkr(1,:,:).*dPcosth0;    %m=0

% ph-component
% GRASP manual signs
% np_mmn_ph_num = -dhkr.*1i.*(-Mmat).*Pcosth.*ephm;   %-m
% np_mpn_ph_num = -dhkr.*1i.*(+Mmat).*Pcosth.*ephp;   %-m
% Hansen/FEKO signs
np_mmn_ph_num = Min1Fact.*dhkr.*1i.*(-Mmat).*Pcosth.*ephm;   %-m
np_mpn_ph_num = Min1Fact.*dhkr.*1i.*(+Mmat).*Pcosth.*ephp;   %-m
% No need for m=0 case - all zeros and handled later
np_mmn_ph = np_mmn_ph_num./sin(th);
np_mpn_ph = np_mpn_ph_num./sin(th);
% Sort out indeterminate forms by L'Hospital
iNaNm = find(abs(np_mmn_ph_num) < eps & abs(sin(th)) < eps);
iNaNp = find(abs(np_mpn_ph_num) < eps & abs(sin(th)) < eps);
% GRASP manual signs
% np_mmn_ph(iNaNm) = -dhkr(iNaNm).*1i.*(-Mmat(iNaNm)).*dPcosth(iNaNm).*ephm(iNaNm)./cos(th(iNaNm));
% np_mpn_ph(iNaNp) = -dhkr(iNaNp).*1i.*(+Mmat(iNaNp)).*dPcosth(iNaNp).*ephp(iNaNp)./cos(th(iNaNp));
% Hansen/FEKO signs
np_mmn_ph(iNaNm) = Min1Fact.*dhkr(iNaNm).*1i.*(-Mmat(iNaNm)).*dPcosth(iNaNm).*ephm(iNaNm)./cos(th(iNaNm));
np_mpn_ph(iNaNp) = Min1Fact.*dhkr(iNaNp).*1i.*(+Mmat(iNaNp)).*dPcosth(iNaNp).*ephp(iNaNp)./cos(th(iNaNp));
clear np_mmn_th_num np_mpn_th_num
clear iNaNm iNaNp

%% Fill the F matrices
Ffact0 = 1/sqrt(2*pi)*1./sqrt(Nmat.*(Nmat+1));
Ffactm = Ffact0.*(+Mmat./Mmat).^(-Mmat);
Ffactp = Ffact0.*(-Mmat./Mmat).^(+Mmat);
% Resize the m=0 factor
Ffact0 = Ffact0(1,:,:);

% Pre-allocate
[Fsmmn,Fsmpn] = deal(zeros(2,MMAX,NMAX,NDIR,3));
[Fsm0n] = deal(zeros(2,1,NMAX,NDIR,3));
% r-components
Fsmmn(1,:,:,:,1) = zeros(1,MMAX,NMAX,NDIR,1);
Fsmmn(2,:,:,:,1) = Ffactm.*np_mmn_r;
Fsmpn(1,:,:,:,1) = zeros(1,MMAX,NMAX,NDIR,1);
Fsmpn(2,:,:,:,1) = Ffactp.*np_mpn_r;
Fsm0n(1,1,:,:,1) = zeros(1,1,NMAX,NDIR,1);
Fsm0n(2,1,:,:,1) = Ffact0.*np_m0n_r;

% th-components
Fsmmn(1,:,:,:,2) = Ffactm.*mp_mmn_th;
Fsmmn(2,:,:,:,2) = Ffactm.*np_mmn_th;
Fsmpn(1,:,:,:,2) = Ffactp.*mp_mpn_th;
Fsmpn(2,:,:,:,2) = Ffactp.*np_mpn_th;
Fsm0n(1,1,:,:,2) = zeros(1,1,NMAX,NDIR,1);
Fsm0n(2,1,:,:,2) = Ffact0.*np_m0n_th;

% ph-components
Fsmmn(1,:,:,:,3) = Ffactm.*mp_mmn_ph;
Fsmmn(2,:,:,:,3) = Ffactm.*np_mmn_ph;
Fsmpn(1,:,:,:,3) = Ffactp.*mp_mpn_ph;
Fsmpn(2,:,:,:,3) = Ffactp.*np_mpn_ph;
Fsm0n(1,1,:,:,3) = Ffact0.*mp_m0n_ph;
Fsm0n(2,1,:,:,3) = zeros(1,1,NMAX,NDIR,1);

% keyboard;



