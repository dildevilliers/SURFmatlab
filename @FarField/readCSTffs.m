function [FF] = readCSTffs(pathName)

% function [FF] = readCSTffs(pathName)
% Loads a CST generated farfield source file pathName.ffs.
%
%
% Inputs:
% path_name - Full path and filename (no extension) string
% Outputs:
% FF - standard farfield object
%
% Dirk de Villiers
% Created: 2015-02-03
% Last edit: 2019-18-02


%Open the data file
% global fid;
[fid, message] = fopen([pathName,'.ffs']);
if (fid==-1)
    error(['Unable to open data file ', fileName, '!']);
end

%===================================================================
% LOAD DATA
%===================================================================

% Read the main header info
freqMarker = '// #Frequencies';
posMarker = '// Position';
zAxisMarker = '// zAxis';
xAxisMarker = '// xAxis';
powerFreqMarker = '// Radiated/Accepted/Stimulated Power , Frequency';
NphNthMarker = '// >> Total #phi samples, total #theta samples';
fieldMarker = '// >> Phi, Theta, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi):';

fCount = 0;
read = 1;
while read
    a = fgetl(fid);
    if strcmp(a,freqMarker) % Read the number of frequencies
        Nf = fscanf(fid,'%i',1);
    end
    if strcmp(a,posMarker) % Read the position
        pos = fscanf(fid,'%f%f%f',[3,1]);
    end
    if strcmp(a,zAxisMarker) % Read the zAxis
        zAxis = fscanf(fid,'%f%f%f',[3,1]);
    end
    if strcmp(a,xAxisMarker) % Read the xAxis
        xAxis = fscanf(fid,'%f%f%f',[3,1]);
    end
    if strncmp(a,powerFreqMarker,length(powerFreqMarker))
        PF = fscanf(fid,'%f',[Nf*4,1]);
        [Prad,Pacc,Pstim,freq] = deal(zeros(1,Nf));
        for ii = 1:Nf
            Prad(ii) = PF((ii-1)*4+1);
            Pacc(ii) = PF((ii-1)*4+2);
            Pstim(ii) = PF((ii-1)*4+3);
            freq(ii) = PF((ii-1)*4+4);
        end
    end
    if strncmp(a,NphNthMarker,length(NphNthMarker))
        NphNth = fscanf(fid,'%i %i',[2,1]);
        Nph = NphNth(1);
        Nth = NphNth(2);
        if fCount == 0
            % Initialise matrices
            [th,ph,Eth,Eph] = deal(zeros(Nth*Nph,Nf));
        end
    end
    
    if strncmp(a,fieldMarker,length(fieldMarker))
        fCount = fCount + 1;
        fDataForm = '%f %f %f %f %f %f';
        fData = fscanf(fid,fDataForm,[6,Nth*Nph])';
        th(:,fCount) = deg2rad(fData(:,2));
        ph(:,fCount) = deg2rad(fData(:,1));
        Eth(:,fCount) = fData(:,3) + 1i.*fData(:,4);
        Eph(:,fCount) = fData(:,5) + 1i.*fData(:,6);
        if fCount == Nf
            read = 0;
        end
    end
end

fclose(fid);

%% Build the FF obj


x = ph(:,1);
y = th(:,1);

E1 = Eth;
E2 = Eph;
E3 = zeros(size(E1));

radEff = Prad./Pacc;
coorType = 'spherical';
polType = 'linear';
gridType = 'PhTh';
freqUnit = 'Hz';


FF = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorType,polType,gridType,freqUnit);

end
