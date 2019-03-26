classdef SWE < FarFieldExpansion
    
    %Class: SWE (Spherical Wave Expansion)
    %Superclass: FarFieldExpansion (abstract)
    %Subclass(es): none
    %Class that represents a set of far-field basis functions from
    %Spherical Wave Expanion principles, based largely on the book
    %"Spherical Near Field Antenna Measurements" (Hansen) and the GRASP
    %technical documentation. Methods allow transformation of a target
    %far-field to a set of coefficients according to the object basis, and
    %vice-versa.
    
    properties
        MMAX(1,1) double {mustBeReal, mustBeFinite}
        NMAX(1,1) double {mustBeReal, mustBeFinite}
        r0(1,1) {mustBeReal}
        freq(:,1) {mustBeReal, mustBeFinite, mustBePositive}
        basis %j-indexed spherical Q-modes represented as FarField objects
        nBasis %number of Q-modes in basis
    end
    
    properties (Constant = true, Hidden = true)
        c0 = physconst('Lightspeed');
        eps0 = 8.854187817000001e-12;
        mu0 = 1.256637061435917e-06;
        eta0 = 3.767303134749689e+02;
    end
    
    methods
        
        function obj = SWE(ang,modes,modestype,freq)
            
            %function obj = SWE(ang,modes,modestype,freq)
            %Constructor method for the SWE object.
            %Inputs:
            %--ang: Input specifying angular range over which basis must be
            %defined (and optionally spherical radius). May either be given
            %as [Nang x 3] vector of [r th ph] values, or as a FarField
            %object.
            %--modes: Mode specifier, either through a minimum sphere
            %radius, a [2 x 1] vector of MMAX/NMAX values or a vector list
            %of Q-mode j-indices.
            %--modestype: Mode type identifier. Can be 'r0' (minimum
            %sphere), 'MNmax' (maximum azimuthal/polar modenumbers) or
            %'list' (j-index list).
            %freq: scalar frequency (in Hz) at which SWE basis should be
            %obtained. Note that SWE objects and their bases are defined at a single
            %frequency only, although they may be used to
            %decompose/reconstruct target far-fields defined over multiple
            %frequencies and frequencies different from the SWE object's frequency.
            %The notion of 'frequency' only relates to the wavenumber value
            %required for parts of the Q-mode calculations (see
            %FsmnFast.m). Users who wish to have basis function sets that
            %vary with frequency (e.g. larger sets of modes as frequency
            %increases for an electrically large antenna) should define one
            %SWE object per frequency that the target far-fields are
            %defined at.
            
            if strcmp(class(ang),'FarField')
                Nang = ang.Nang;
                r = inf.*ones(size(ang.th));
                th = ang.th;
                ph = ang.ph;
                if nargin < 4 %Default to highest frequency in FFobj if frequency is not explicitly specified
                    freq = ang.freq(end);
                end
            else
                Nang = size(ang,1);
                r = ang(:,1);
                th = ang(:,2);
                ph = ang(:,3);
                if nargin < 4
                    disp('No frequency specified - defaulting to 1 GHz')
                    freq = 1e9;
                end
            end
            
            %handle maximum numbers of modes
            switch modestype
                case 'r0'
                    r0 = modes;
                    lam = SWE.c0/freq;
                    k = 2*pi/lam;
                    MMAXin = ceil(k*r0 + max(3.6*(k*r0)^(1/3),10));
                    %     MMAXin = 2.*floor((k.*r0 + 3.*((k.*r0)^(1./3)))./2);
                    NMAXin = MMAXin;
                case 'MNmax'
                    r0 = 0;
                    MMAXin = modes(1);
                    NMAXin = modes(2);
                case 'list'
                    error('j-index list option not yet implemented')
            end
            
            [Fsm0n,Fsmmn,Fsmpn] = FsmnFast(r,th,ph,MMAXin,NMAXin,freq);
            
            %Pack each basis mode into its own FarField object
            basis = cell(2*NMAXin*(NMAXin+2),1);
            for jcurr = 1:2*NMAXin*(NMAXin+2)
                smn = j2smn(jcurr);
                [scurr mcurr ncurr] = deal(smn(1),smn(2),smn(3));
                if abs(mcurr) <= MMAXin && ncurr <= NMAXin %only fill in basis if its smn-index falls within the max mode numbers
                    if mcurr == 0
                        Fr(:,1) = reshape(Fsm0n(scurr,1,ncurr,:,1),Nang,1);
                        Fth(:,1) = reshape(Fsm0n(scurr,1,ncurr,:,2),Nang,1);
                        Fph(:,1) = reshape(Fsm0n(scurr,1,ncurr,:,3),Nang,1);
                    elseif mcurr < 0
                        Fr(:,1) = reshape(Fsmmn(scurr,-mcurr,ncurr,:,1),Nang,1);
                        Fth(:,1) = reshape(Fsmmn(scurr,-mcurr,ncurr,:,2),Nang,1);
                        Fph(:,1) = reshape(Fsmmn(scurr,-mcurr,ncurr,:,3),Nang,1);
                    else
                        Fr(:,1) = reshape(Fsmpn(scurr,mcurr,ncurr,:,1),Nang,1);
                        Fth(:,1) = reshape(Fsmpn(scurr,mcurr,ncurr,:,2),Nang,1);
                        Fph(:,1) = reshape(Fsmpn(scurr,mcurr,ncurr,:,3),Nang,1);
                    end
                    basis{jcurr} = FarField(ph,th,Fth,Fph,Fr,freq);
                end
            end
            
            %Set all properties for object from SWE operation
            obj.NMAX = NMAXin;
            obj.MMAX = MMAXin;
            obj.r0 = r0;
            obj.basis = basis;
            obj.nBasis = length(basis);
            obj.freq = freq;
            
            %Validate properties inherited from abstract class (can't be validated inline like other properties...)
            mustBeReal(obj.nBasis);
            mustBeFinite(obj.nBasis);
            mustBePositive(obj.nBasis);
            assert(iscell(obj.basis),'Basis property must be cell array');
            assert(all(size(obj.basis) == [1 obj.nBasis]) || all(size(obj.basis) == [obj.nBasis 1]),'Basis property must be [nBasis x 1] or [1 x nBasis] 1D cell array');
            for bb = 1:numel(obj.basis)
                assert(isa(obj.basis{bb},'FarField') || isempty(obj.basis{bb}),'All Basis property elements must be FarField objects')
            end
            
        end
    end
    
    methods (Static)
        
        function [Qjout,Pout] = farField2Expansion(SWEobj,FFobj,iBasis)
            
            %function [Qjout,Pout] = farField2Expansion(SWEobj,FFobj,iBasis)
            %SWE concrete implementation of expanding a target FarField
            %object into SWE Q-coefficients and radiated power terms.
            %Inputs:
            %--SWEobj: SWE object to perform expansion with.
            %--FFobj: target FarField object to expand.
            %--iBasis: optional vector of basis indices to use in the
            %expansion - all basis functions not selected by their indices
            %are masked out (zero'd) in the output expansion.
            
            %Get this SWE object's basis functions
            if nargin < 3
                currBasis = SWEobj.basis;
            else
                currBasis = SWEobj.basis{iBasis};
            end
            
            %extract farfield vectors from FFobj, build farfield struct for FF2SWE function
            [Eth, Eph, ~] = FFobj.getEspherical;
            FFstr.Eth = Eth;
            FFstr.Eph = Eph;
            FFstr.th = FFobj.th;
            FFstr.ph = FFobj.ph;
            FFstr.Nth = length(unique(FFobj.th));
            FFstr.Nph = length(unique(FFobj.ph));
            FFstr.freq = FFobj.freq;
            FFstr.Nf = FFobj.Nf;
            
            %Turn basis FarField objects into a useable F struct for FF2SWE
            Fsm0n = zeros(2,1,SWEobj.NMAX,SWEobj.basis{1}.Nang,3);
            Fsmmn = zeros(2,SWEobj.MMAX,SWEobj.NMAX,SWEobj.basis{1}.Nang,3);
            Fsmpn = zeros(2,SWEobj.MMAX,SWEobj.NMAX,SWEobj.basis{1}.Nang,3);
            for jj = 1:SWEobj.nBasis
                smn = j2smn(jj);
                [scurr mcurr ncurr] = deal(smn(1),smn(2),smn(3));
                if abs(mcurr) <= abs(SWEobj.MMAX) && ncurr <= SWEobj.NMAX
                    if mcurr == 0
                        Fsm0n(scurr,1,ncurr,:,1) = SWEobj.basis{jj}.E3;
                        Fsm0n(scurr,1,ncurr,:,2) = SWEobj.basis{jj}.E1;
                        Fsm0n(scurr,1,ncurr,:,3) = SWEobj.basis{jj}.E2;
                    elseif mcurr < 0
                        Fsmmn(scurr,-mcurr,ncurr,:,1) = SWEobj.basis{jj}.E3;
                        Fsmmn(scurr,-mcurr,ncurr,:,2) = SWEobj.basis{jj}.E1;
                        Fsmmn(scurr,-mcurr,ncurr,:,3) = SWEobj.basis{jj}.E2;
                    else
                        Fsmpn(scurr,mcurr,ncurr,:,1) = SWEobj.basis{jj}.E3;
                        Fsmpn(scurr,mcurr,ncurr,:,2) = SWEobj.basis{jj}.E1;
                        Fsmpn(scurr,mcurr,ncurr,:,3) = SWEobj.basis{jj}.E2;
                    end
                end
            end
            Fin = struct('Fsm0n',Fsm0n,'Fsmmn',Fsmmn,'Fsmpn',Fsmpn);
            Fincell{1} = Fin;
            %perform SWE on farfield vectors, get Q-modes and basis
            [~,Qjout,Pout,~] = FF2SWE(FFstr,SWEobj.NMAX,SWEobj.MMAX,[],0,Fincell);
        end
        
        function FFobj = expansion2FarField(SWEobj,Qj,P,iBasis)
            
            %function FFobj = expansion2FarField(SWEobj,Qj,P,iBasis)
            %SWE concrete implementation of reconstructing a far-field from
            %SWE expansion coefficients and powers.
            %Inputs:
            %--SWEobj: SWE object to perform reconstruction with
            %--Qj: Q-coefficients struct, in the format given by farField2Expansion. 
            %--P: optional vector of power terms, in the format given by farField2Expansion. 
            %--iBasis: optional vector of basis indices to use in the
            %reconstruction - all basis functions not selected by their indices
            %are masked out (zero'd) in the reconstructed output FarField object.
            
            %Set spherical grid values
            th = SWEobj.basis{1}.th;
            ph = SWEobj.basis{1}.ph;
            r = inf.*ones(size(th));
            
            %If a set of basis indices were given as 4th input, use them to mask Q-coefficients
            Qjin = Qj;
            if nargin > 3
                Qmask = zeros(SWEobj.nBasis,1);
                Qmask(iBasis) = 1;
                Qjin.Q = Qjin.Q.*Qmask;
            end
            
            %Pack basis functions into matrix for SWE2FF
            F = zeros(2*SWEobj.basis{1}.Nang,SWEobj.nBasis);
            for bb = 1:SWEobj.nBasis
                if ~isempty(SWEobj.basis{bb})
                    F(:,bb) = [SWEobj.basis{bb}.E1;SWEobj.basis{bb}.E2];
                end
            end
            
            %Reconstruct pattern and pack it into FarField object output
            if nargin < 3
                FFstr = SWE2FF(Qjin,r,th,ph,F);
                FFobj = FarField(FFstr.ph,FFstr.th,FFstr.Eth,FFstr.Eph,zeros(size(FFstr.Eth)),FFstr.freq);
            else
                FFstr = SWE2FF(Qjin,r,th,ph,F,[P.P]);
                FFobj = FarField(FFstr.ph,FFstr.th,FFstr.Eth,FFstr.Eph,zeros(size(FFstr.Eth)),FFstr.freq,FFstr.Prad);
            end
            
        end
        
        
    end
end