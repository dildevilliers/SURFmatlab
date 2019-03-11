classdef SWE < FarFieldExpansion
    
    properties
        NMAX(1,1) double {mustBeReal, mustBeFinite}
        MMAX(1,1) double {mustBeReal, mustBeFinite}
        r0(1,1) {mustBeReal}
        freq(:,1) {mustBeReal, mustBeFinite, mustBePositive}
        basis %j-indexed spherical Q-modes
        nBasis
    end
    
    properties (Constant = true, Hidden = true)
        c0 = physconst('Lightspeed');
        eps0 = 8.854187817000001e-12;
        mu0 = 1.256637061435917e-06;
        eta0 = 3.767303134749689e+02;
    end
    
    methods
        
        function obj = SWE(FFobj,modes,modestype)
            
            %handle maximum numbers of modes
            switch modestype
                case 'r0'
                    assert(length(modes) == FFobj.Nf, 'number of r0 values must equal FarField object Nf');
                    r0 = modes;
                    lam = SWE.c0./FFobj.freq;
                    k = 2.*pi./lam;
                    MMAXin = ceil(k.*r0 + max(3.6.*(k.*r0).^(1/3),10));
                    %                     MMAXin = 2.*floor((k.*r0 + 3.*((k.*r0)^(1./3)))./2);
                    NMAXin = MMAXin;
                case 'MNmax'
                    assert(size(modes,1) == FFobj.Nf, 'number of MMAX/NMAX values must equal FarField object Nf');
                    r0 = zeros(size(modes),1);
                    MMAXin = modes(:,1);
                    NMAXin = modes(:,2);
            end
            
            for ff = 1:1%FFobj.Nf
                [Fsm0n{ff},Fsmmn{ff},Fsmpn{ff}] = FsmnFast(inf.*ones(size(FFobj.th)),FFobj.th,FFobj.ph,MMAXin(ff),NMAXin(ff),FFobj.freq(ff));
            end
            
            %Pack each basis mode into its own FarField object
            for jcurr = 1:2*size(Fsmmn{ff},3)*(size(Fsmmn{ff},3)+2)
                smn = j2smn(jcurr);
                [scurr mcurr ncurr] = deal(smn(1),smn(2),smn(3));
                
                for ff = 1:1%FFobj.Nf
                    if mcurr == 0
                        Fr(:,ff) = reshape(Fsm0n{ff}(scurr,1,ncurr,:,1),FFobj.Nang,1);
                        Fth(:,ff) = reshape(Fsm0n{ff}(scurr,1,ncurr,:,2),FFobj.Nang,1);
                        Fph(:,ff) = reshape(Fsm0n{ff}(scurr,1,ncurr,:,3),FFobj.Nang,1);
                    elseif mcurr < 0
                        Fr(:,ff) = reshape(Fsmmn{ff}(scurr,-mcurr,ncurr,:,1),FFobj.Nang,1);
                        Fth(:,ff) = reshape(Fsmmn{ff}(scurr,-mcurr,ncurr,:,2),FFobj.Nang,1);
                        Fph(:,ff) = reshape(Fsmmn{ff}(scurr,-mcurr,ncurr,:,3),FFobj.Nang,1);
                    else
                        Fr(:,ff) = reshape(Fsmpn{ff}(scurr,mcurr,ncurr,:,1),FFobj.Nang,1);
                        Fth(:,ff) = reshape(Fsmpn{ff}(scurr,mcurr,ncurr,:,2),FFobj.Nang,1);
                        Fph(:,ff) = reshape(Fsmpn{ff}(scurr,mcurr,ncurr,:,3),FFobj.Nang,1);
                    end
                end
               
                basis(jcurr) = FarField(FFobj.ph,FFobj.th,Fth,Fph,Fr,FFobj.freq);
            end
            
            
            %Set all properties for object from SWE operation
            obj.NMAX = NMAXin;
            obj.MMAX = MMAXin;
            obj.r0 = r0;
            obj.basis = basis;
            obj.nBasis = length(basis);
            obj.freq = FFobj.freq;
        end
        
        
        %         function obj = plotCoeffs(obj,iCoeffs)
        %
        %             if (nargin < 2)
        %
        %             end
        %
        %             p = inputparser;
        %             addParameter(p,'outputType','realimag');
        %             addParameter(p,'scaleMag','lin');
        %
        %
        %             figure
        %             switch plotType
        %                 case ('mag')
        %                     plot(abs(obj.coeffs(iCoeffs)));
        %                 case('phase')
        %                     plot(deg2rad(angle(obj.coeffs(iCoeffs))));
        %                 case('real')
        %                     plot(real(obj.coeffs(iCoeffs)));
        %                 case('imag')
        %                     plot(imag(obj.coeffs(iCoeffs)));
        %                 case('realimag')
        %                     subplot(1,2,1)
        %                     plot(real(obj.coeffs(iCoeffs)));
        %                     subplot(1,2,1)
        %                     plot(imag(obj.coeffs(iCoeffs)));
        %                 otherwise
        %                     error('plotType argument unrecognised')
        %             end
        %         end
        
        
    end
    
    methods (Static)
        
        function [Qjout,Pout,Fout] = farField2Expansion(obj,FFobj,iBasis)
            
            %Get this SWE object's basis functions
            if nargin < 3
                currBasis = obj.basis;
            else
                currBasis = obj.basis(iBasis);
            end
            
            for ii = 1:length(currBasis)
                Fin = obj.basis;
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
            
            %perform SWE on farfield vectors, get Q-modes and basis
            [~,Qjout,Pout,Fout] = FF2SWE(FFstr,obj.NMAX,obj.MMAX,[],0);
            
        end
        
        function FFobj = expansion2FarField(obj,Qj,P)
            
            th = obj.basis(1).th;
            ph = obj.basis(1).ph;
            r = inf.*ones(size(th));
            for bb = 1:length(obj.basis)
                F(:,bb) = [obj.basis(bb).E1;obj.basis(bb).E2;obj.basis(bb).E3];
            end
            %Get weights and basis functions
            %Resconstruct
            if nargin < 3
                FFstr = SWE2FF(Qj,r,th,ph,F);
                FFobj = FarField(FFstr.ph,FFstr.th,FFstr.Eth,FFstr.Eph,zeros(size(FFstr.Eth)),FFstr.freq);
            else
                FFstr = SWE2FF(Qj,r,th,ph,F,P.P);
                FFobj = FarField(FFstr.ph,FFstr.th,FFstr.Eth,FFstr.Eph,zeros(size(FFstr.Eth)),FFstr.freq,FFstr.Prad);
            end

        end
        
    end
    
    
    
end