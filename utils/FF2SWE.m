function [Qsmnout,Qjout,Pout,Fout] = FF2SWE(FFin,NMAXin,MMAXin,opts,r0,Fin)
% function [Qsmnout,Pout,Qjout,Fout] = FF2SWE(FFin,NMAXin,MMAXin,r0,opts)
% Performs Spherical Wave Expansion (SWE) on a set of input far/nearfields (defined
% in spherical coordinates), and returns the Q-expansion coefficients and
% related power information.
% Inputs:
%   --FFin - NX-by-1 struct array of FF structures of length Nx, each element
%   containing farfield evaluated at NF frequencies
%   --NMAXin - User-defined maximum number of polar modes to use during SWE.
%   Can be given as scalar (applied uniformly across all FFin's elements and frequencies),
%   NX-by-1 vector (applied element-wise to FFin, 1-by-NF vector (applied uniformly across all FFin's elements but
%   element-wise in frequency), NX-by-NF array or empty (calculates maximum
%   modes numbers using the minimum sphere).
%   --MMAXin - User-defined maximum number of azimuthal modes to use during
%   SWE. Can be given in the same ways as NMAXin.
%   --r0 - Optional radii of minimum spheres enclosing each antenna associated with
%   the fields given by FFin, in metres. Can be given as a scalar (assumes
%   all antennas associated with FFin are enclosed by the single given
%   sphere radius), or as an NX-by-1 vector. If NMAXin/MMAXin are not fully
%   specified for all FFin/frequencies, this must be set
%   --opts - optional options for the method:
%       opts.fekor0 - if 1, uses (what appears to be) FEKO's method of
%       obtaining NMAX
% Returns:
%   --Qsmnout - NX-by-NF struct array, each element containing fields as that
%   of the output to readSPH.m
%   --Qjout - As Qsmnout, except that instead of the 3 smn-indexed
%   Q-coefficient fields, a single field Q is given that provides the
%   Q-coefficients with j-indexing (as described in Hansen)
%   --Pout - NX-by-NF struct array containing POWERM (as from readSPH.m
%   output) as well as the total radiated power P
%   --Fout - NX-by-NF struct array of the SWE modal basis functions used to
%   obtain the Q-coefficients

% Created: 2017-04-19
% Updated: 2017-10-26
% Dirk de Villiers, Brandt Klopper


%% Constants

c0 = 299792458;
eta0 = 3.767303134749689e+02;

%% Basic error checking

% if (FFin(1).Nf > 1) assert(isempty(FFin.freq) || all(ismember(reshape([FFin.freq],length(FFin(1).freq),length(FFin))',FFin(1).freq,'rows')),'All FF.freq fields should be identical'); end
if (FFin(1).Nf > 1) assert(all(ismember(reshape([FFin.freq],length(FFin(1).freq),length(FFin))',FFin(1).freq,'rows')),'All FF.freq fields should be identical'); end
assert(all([FFin.Nth] == FFin(1).Nth),'All FF.Nth fields should be identical')
assert(all([FFin.Nph] == FFin(1).Nph),'All FF.Nph fields should be identical')
if (nargin > 4) assert(length(r0) == length(FFin) || length(r0) == 1,'r0 should be scalar or a vector of the same length as FFin'); end

%% Define spherical coordinate variable ranges

th = FFin(1).th(:,1);
ph = FFin(1).ph(:,1);
if isfield(FFin,'r')
    r = FFin(1).r(:,1);
else
    r = ones(size(th)).*inf;
end

%% Other Preliminaries

NX = length(FFin);
[NDIR,NF] = size(FFin(1).Eth);

if isfield(FFin,'freq') && ~isempty(FFin(1).freq)
    freq = reshape(FFin(1).freq,1,NF);
    powerNorm = false;
else
    warning('No frequency found - normalizing power in modes to 0.5 W');
    freq = ones(1,NF);
    powerNorm = true;
end
lam = c0./freq;
k = 2.*pi./lam;

if isfield(FFin,'Er')
    getR = 1;
else
    getR = 0;
end

%Handle optional input r0
if nargin < 4
    r0 = 1;
    opts = [];
elseif nargin < 5
    r0 = 1;
end
r0  = r0.*ones(NX,1); %In the case that r0 is a scalar, this will simply make it a vector that is compatible with the subsequent loops

%% Fill Matrix

getallF = 1; %flag to signify that F-matrices should be calculated for each iteration aa as well as ff along frequency axis of each element of FFin

[Qsm0n,Qsmmn,Qsmpn,POWERM,P] = deal(cell(NX,NF));

if nargin < 6
    if getR == 0 && length(NMAXin) == 1 && length(MMAXin) == 1 %Farfield case where NMAX, MMAX are uniform across input set
        [F0all,Fmall,Fpall] = FsmnFast(r,th,ph,MMAXin,NMAXin,1); %frequency set to 1 because farfield F-functions are independent of f
        [Fsm0n{1:NX,1:NF}] = deal(F0all);
        [Fsmmn{1:NX,1:NF}] = deal(Fmall);
        [Fsmpn{1:NX,1:NF}] = deal(Fpall);
        getallF = 0;
        disp('Uniform mode limits detected- only calculating one F-matrix for all inputs')
    end
else
    if getR == 0 && length(NMAXin) == 1 && length(MMAXin) == 1
        [Fsm0n{1:NX,1:NF}] = deal(Fin{1}.Fsm0n);
        [Fsmmn{1:NX,1:NF}] = deal(Fin{1}.Fsmmn);
        [Fsmpn{1:NX,1:NF}] = deal(Fin{1}.Fsmpn);
    else
        for xx = 1:NX
            for ff = 1:NF
                Fsm0n{xx,ff} = Fin{xx,ff}.Fsm0n;
                Fsmmn{xx,ff} = Fin{xx,ff}.Fsmmn;
                Fsmpn{xx,ff} = Fin{xx,ff}.Fsmpn;
            end
        end
    end
end


for xx = 1:NX
    for ff = 1:NF
        if isfield(opts,'fekor0') && opts.fekor0 == 1
            N =  2*floor((k(ff)*r0(xx) + 3*((k(ff)*r0(xx))^(1/3)))/2); % EXPERIMENTAL - based on FEKO docs and patterns observed in .sph files from FEKO
        else
            N = ceil(k(ff)*r0(xx) + max(3.6*(k(ff)*r0(xx)).^(1/3),10));
        end
        if isempty(NMAXin)
            NMAX{xx,ff} = N;
        elseif length(NMAXin) == 1
            NMAX{xx,ff} = NMAXin;
        elseif size(NMAXin,1) == 1
            NMAX{xx,ff} = NMAXin(ff);  % Will throw error for wrong lengths
        elseif size(NMAXin,2) == 1
            NMAX{xx,ff} = NMAXin(xx);  % Will throw error for wrong lengths
        else
            NMAX{xx,ff} = NMAXin(xx,ff);  % Will throw error for wrong lengths
        end
        if isempty(MMAXin)
            MMAX{xx,ff} = N;
        elseif length(MMAXin) == 1
            MMAX{xx,ff} = MMAXin;
        elseif size(MMAXin,1) == 1
            MMAX{xx,ff} = MMAXin(ff);  % Will throw error for wrong lengths
        elseif size(MMAXin,2) == 1
            MMAX{xx,ff} = MMAXin(xx);  % Will throw error for wrong lengths
        else
            MMAX{xx,ff} = MMAXin(xx,ff);  % Will throw error for wrong lengths
        end
        
        if getallF == 1 && nargin < 6  % Nearfield case, or farfields with variations in NMAXin, MMAxin along xx and ff
            [Fsm0n{xx,ff},Fsmmn{xx,ff},Fsmpn{xx,ff}] = FsmnFast(r,th,ph,MMAX{xx,ff},NMAX{xx,ff},freq(ff));
        end
        
        % Build the matrix for different components first
        % First get indices of valid mode numbers
        [SS,MM,NN] = ndgrid(1:2,1:MMAX{xx,ff},1:NMAX{xx,ff});
        NN(NN<MM) = 0;
        validModePos = find(NN>0);
        NvalidMode = length(validModePos);
        [Fr,Fth,Fph] = deal(zeros(NDIR,2*(NMAX{xx,ff}+NvalidMode)));
        for dd = 1:NDIR
            if getR == 1
                Fr0 = reshape(Fsm0n{xx,ff}(:,:,:,dd,1),2*NMAX{xx,ff},1);
                Frm = reshape(Fsmmn{xx,ff}(:,:,:,dd,1),2*MMAX{xx,ff}*NMAX{xx,ff},1);
                Frp = reshape(Fsmpn{xx,ff}(:,:,:,dd,1),2*MMAX{xx,ff}*NMAX{xx,ff},1);
                Fr(dd,:) = [Fr0;Frm(validModePos);Frp(validModePos)].';
            else
                Fr = [];
            end
            
            Fth0 = reshape(Fsm0n{xx,ff}(:,:,:,dd,2),2*NMAX{xx,ff},1);
            Fthm = reshape(Fsmmn{xx,ff}(:,:,:,dd,2),2*MMAX{xx,ff}*NMAX{xx,ff},1);
            Fthp = reshape(Fsmpn{xx,ff}(:,:,:,dd,2),2*MMAX{xx,ff}*NMAX{xx,ff},1);
            Fth(dd,:) = [Fth0;Fthm(validModePos);Fthp(validModePos)].';

            Fph0 = reshape(Fsm0n{xx,ff}(:,:,:,dd,3),2*NMAX{xx,ff},1);
            Fphm = reshape(Fsmmn{xx,ff}(:,:,:,dd,3),2*MMAX{xx,ff}*NMAX{xx,ff},1);
            Fphp = reshape(Fsmpn{xx,ff}(:,:,:,dd,3),2*MMAX{xx,ff}*NMAX{xx,ff},1);
            Fph(dd,:) = [Fph0;Fphm(validModePos);Fphp(validModePos)].';
            
        end
        F = [Fr;Fth;Fph];   % This is only populated for the valid mode columns - wont reshape to rectangular matrices
        
        % RHS vector
        if getR
            E = [FFin(xx).Er(:,ff);FFin(xx).Eth(:,ff);FFin(xx).Eph(:,ff)];
        else
            E = [FFin(xx).Eth(:,ff);FFin(xx).Eph(:,ff)];
        end
        
        % Solve for Q - NB: Qv is divided by sqrt(eta0) to produce Q-coefficients with same magnitude as those provided by FEKO/GRASP
        keyboard
        Qv = (F\E)./sqrt(eta0);   % This is for all the valid modes (not for m > n modes which are included in the matrices from FsmnFast.m)

        
        % Expand and repack Q into the standard matrices used by SWEgetField.m
        % zero m
        Qm0 = Qv(1:2*NMAX{xx,ff});
        Qsm0n{xx,ff} = reshape(Qm0,2,1,NMAX{xx,ff});
        [Qmm,Qmp] = deal(zeros(2*MMAX{xx,ff}*NMAX{xx,ff},1));
        % negative m
        Qmm(validModePos) = Qv(2*NMAX{xx,ff}+1:2*NMAX{xx,ff}+NvalidMode);
        Qsmmn{xx,ff} = reshape(Qmm,2,MMAX{xx,ff},NMAX{xx,ff});
        % positive m
        Qmp(validModePos) = Qv(2*NMAX{xx,ff}+NvalidMode+1:end);
        Qsmpn{xx,ff} = reshape(Qmp,2,MMAX{xx,ff},NMAX{xx,ff});
        
        % Calculate the power in each of the m-modes for normalization
        POWERM{xx,ff}(1,1) = 0;
        POWERM{xx,ff}(1,2) = sum(sum(abs(Qsm0n{xx,ff}).^2,3),1);
        POWERM{xx,ff}(2:MMAX{xx,ff}+1,1) = 1:MMAX{xx,ff};
        POWERM{xx,ff}(2:MMAX{xx,ff}+1,2) = sum(sum(abs(Qsmmn{xx,ff}).^2,3),1) + sum(sum(abs(Qsmpn{xx,ff}).^2,3),1);
        
        if powerNorm
            PnormFact = 0.5/max(POWERM{xx,ff}(:,2));
            POWERM{xx,ff}(:,2) = POWERM{xx,ff}(:,2)*PnormFact;
            Qsm0n{xx,ff} = Qsm0n{xx,ff}.*sqrt(PnormFact);
            Qsmmn{xx,ff} = Qsmmn{xx,ff}.*sqrt(PnormFact);
            Qsmpn{xx,ff} = Qsmpn{xx,ff}.*sqrt(PnormFact);
        end
        
        %Calculate total radiated power via sum of squares of the modal coefficients
        P{xx,ff} = 0.5*(sum(sum(sum(abs(Qsm0n{xx,ff}).^2,3),2)) + sum(sum(sum(abs(Qsmmn{xx,ff}).^2,3),2)) + sum(sum(sum(abs(Qsmpn{xx,ff}).^2,3),2)));
        
    end
end

%% Output structures

Qsmnout = struct('Qsm0n',Qsm0n,'Qsmmn',Qsmmn,'Qsmpn',Qsmpn,'freq',repmat(num2cell(freq),NX,1),'MMAX',MMAX,'NMAX',NMAX);
Qjout = Qsmn2j(Qsmnout);
Fout = struct('Fsm0n',Fsm0n,'Fsmmn',Fsmmn,'Fsmpn',Fsmpn);

if powerNorm
    Pout = struct('POWERM',POWERM,'P',P,'PnormFact',Pnormfact);
else
    Pout = struct('POWERM',POWERM,'P',P);
end
end
