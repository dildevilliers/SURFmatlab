function [Qsmn,Qsmnprime,POWERM,P,Qj,Qjprime] = readSPH(pathName)

% function [Q,Qprime,POWERM,P] = readSPH(pathName)
% Function to read a spherical mode expansion file with extension .sph
% Returns the Q, Qprime and POWERM structures as defined by the GRAPS
% technical description/manual.
% It also returns the total radiated power (sum of modal powers) in P
% Each Q is a structure with:
% Qsm0n: Zero m-modes
% Qsmmn: Negative m-modes
% Qsmnn: Positive m-modes
% freq: 1xNf array of frequency values (could be empty if none are specified)
% MMAX, NMAX: Both functions of frequency containing the number of modes
% for each frequency
% Each of the Q* above are cell arrays of length frequency, with each cell a
% matrix of size (2 x MMAX x NMAX), with MMAX and NMAX generally also a
% function of frequency

[fid, message] = fopen([pathName,'.sph']);

read = true;
fcount = 1; % Frequency counter
lcount = 1; % Line counter - resets after each frequency
while read
%     if lcount == 1
%         disp(['Reading file: ',pathName,'.sph; Frequency point ', num2str(fcount),'...'])
%     end
    a = fgetl(fid);
    if a == -1, break; end
    switch lcount
        case 1
            PRGTAG = a;
        case 2
            IDSTRG = a;
        case 3
            N = sscanf(a,'%i',4);
            NTHE(fcount) = N(1);
            NPHI(fcount) = N(2);
            NMAX(fcount) = N(3);
            MMAX(fcount) = N(4);
        case 4
            TEXT_4 = a;
            aSplit = strsplit(a);
            if strcmp(aSplit{2},'Frequency')
                freq(fcount) = str2double(aSplit{4});
            else
                freq = [];
            end
        case 5
            angles1 = sscanf(a,'%f',5);
        case 6
            angles2 = sscanf(a,'%f',5);
        case 7
            TEXT_DUMMY1 = a;
        case 8
            TEXT_DUMMY2 = a;
            % Initialise reading structure
            Ncoeff(fcount) = NMAX(fcount);
            for mm = 1:MMAX(fcount)
                Ncoeff = Ncoeff + 2*(NMAX(fcount) - (mm - 1));
            end
            Nlines_F(fcount) = Ncoeff(fcount) + MMAX(fcount) + 9;
            SWE_TEXT_1F = cell(Ncoeff(fcount) + MMAX(fcount)+1,1);
        otherwise
            SWE_TEXT_1F{lcount - 8} = a;
            SWE_TEXT{fcount} = SWE_TEXT_1F;
            if lcount >= Nlines_F(fcount)
                fcount = fcount + 1;
                lcount = 0;
            end
    end
    lcount = lcount + 1;
end
fclose(fid);
NF = fcount-1;
POWERM = cell(NF,1);
[Qsm0n_prime,Qsmmn_prime,Qsmpn_prime] = deal(cell(NF,1));
[Qsm0n,Qsmmn,Qsmpn] = deal(cell(NF,1));
P = zeros(1,NF);
% Loop over frequency
for ff = 1:NF
%     disp(['Building coefficient matrix set ', num2str(ff),' of ', num2str(NF), '...'])
    % Initialise mode power matrix (M,POWERM)
    POWERM{ff} = zeros(MMAX(ff)+1,2);
    
    % Store Q-mode blocks
    % m = 0 special case
    % Read the power in the modes
    mStartLine = 1;
    POWERM{ff}(1,:) = sscanf(SWE_TEXT{ff}{mStartLine},'%i%f',2);
    % Read the mode coefficients
    Qsm0n_prime{ff} = zeros(2,1,NMAX(ff));   % m = 0, +m and -m
    for nn = 1:NMAX(ff)
        Q = sscanf(SWE_TEXT{ff}{1+nn},'%f',4);
        Qsm0n_prime{ff}(1,1,nn) = Q(1) + 1i.*Q(2);
        Qsm0n_prime{ff}(2,1,nn) = Q(3) + 1i.*Q(4);
    end
    % Rest of the modes
    [Qsmmn_prime{ff},Qsmpn_prime{ff}] = deal(zeros(2,MMAX(ff),NMAX(ff)));
    mStartLine = mStartLine + NMAX(ff) + 1;
    for mm = 1:MMAX(ff)
        POWERM{ff}(1+mm,:) = sscanf(SWE_TEXT{ff}{mStartLine},'%i%f',2);
        for nn = mm:NMAX(ff)
            Qm_min = sscanf(SWE_TEXT{ff}{mStartLine+2*((nn-mm)+1)-1},'%f',4);
            Qm_plus = sscanf(SWE_TEXT{ff}{mStartLine+2*((nn-mm)+1)},'%f',4);
            Qsmmn_prime{ff}(1,mm,nn) = Qm_min(1) + 1i*Qm_min(2);
            Qsmmn_prime{ff}(2,mm,nn) = Qm_min(3) + 1i*Qm_min(4);
            Qsmpn_prime{ff}(1,mm,nn) = Qm_plus(1) + 1i*Qm_plus(2);
            Qsmpn_prime{ff}(2,mm,nn) = Qm_plus(3) + 1i*Qm_plus(4);
        end
        mStartLine = mStartLine + 2*(NMAX(ff)-(mm-1)) + 1;
    end
    % Sort out normalisation
    Qsm0n{ff} = sqrt(8*pi).*conj(Qsm0n_prime{ff});
    Qsmmn{ff} = sqrt(8*pi).*conj(Qsmmn_prime{ff});
    Qsmpn{ff} = sqrt(8*pi).*conj(Qsmpn_prime{ff});
    % Get total power
    P(ff) = sum(POWERM{ff}(:,2));
    
%     % Test the total power - should be 4*pi for the dBi normalised GRASP cases
%     Pn(ff) = 0.5*(sum(sum(sum(abs(Qsm0n{ff}).^2,3),2)) + sum(sum(sum(abs(Qsmmn{ff}).^2,3),2)) + sum(sum(sum(abs(Qsmpn{ff}).^2,3),2)));
%     % % Hence, below should be 1
%     SumQ(ff) = (sum(sum(sum(abs(Qsm0n_prime{ff}).^2,3),2)) + sum(sum(sum(abs(Qsmmn_prime{ff}).^2,3),2)) + sum(sum(sum(abs(Qsmpn_prime{ff}).^2,3),2)));
    
end

Qsmn = struct('Qsm0n',Qsm0n,'Qsmmn',Qsmmn,'Qsmpn',Qsmpn,'freq',num2cell(freq)','MMAX',num2cell(MMAX)','NMAX',num2cell(NMAX)');
Qsmnprime = struct('Qsm0n',Qsm0n_prime,'Qsmmn',Qsmmn_prime,'Qsmpn',Qsmpn_prime,'freq',num2cell(freq)','MMAX',num2cell(MMAX)','NMAX',num2cell(NMAX)');

Qj = Qsmn2j(Qsmn);
Qjprime = Qsmn2j(Qsmnprime);
%  keyboard

    
