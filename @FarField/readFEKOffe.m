function FF = readFEKOffe(pathName)

%Name: readFEKOffe.m
%Description:
%   Function to create a Farfield object from a FEKO .ffe far-field output
%   file. 
%Inputs:
% th: column vector [Nang x 1] of th angles in rad
% ph: column vector [Nang x 1] of ph angles in rad
%Outputs:
% --FF: Farfield object containing parameters as derived from target .ffe
% file.


%Physical constants
c0 = physconst('Lightspeed');
eps0 = 8.854187817000001e-12;
mu0 = 1.256637061435917e-06;
eta0 = 3.767303134749689e+02;

%Open the data file
[fid, message] = fopen([pathName,'.ffe']);
if (fid==-1)
    error(['Unable to open data file ' fileName '!']);
end

% Read the main header info
coorSysMarker = '#Coordinate System:';
freqMarker = '#Frequency:';
NthMarker = '#No. of Theta Samples:';
NphMarker = '#No. of Phi Samples:';
fieldMarker = '#"Theta""Phi"';


%===================================================================
% LOAD DATA
%===================================================================

fCount = 0;
read = 1;
while read
    a = fgetl(fid);
    if a == -1
        read = 0;
        break;
    end
    if strncmp(a,coorSysMarker,length(coorSysMarker)) % Read the coordinate system type
        coorSysCell = textscan(a,'%s%s%s');
        coorSys = coorSysCell{3};
    end
    if strncmp(a,freqMarker,length(freqMarker)) % Read the number of frequencies
        freqCell = textscan(a,'%s%n');
        freq(fCount+1) = freqCell{2};
    end
    if strncmp(a,NthMarker,length(NthMarker))
        NthCell = textscan(a,'%s%s%s%s%n');
        Nth = NthCell{5};
    end
    if strncmp(a,NphMarker,length(NphMarker))
        NphCell = textscan(a,'%s%s%s%s%n');
        Nph = NphCell{5};
        
        
    end
    
    %     if strncmp(a,fieldMarker,length(fieldMarker))
    aNoSpace = a;
    aNoSpace(ismember(a,' ')) = [];
    if strncmp(aNoSpace,fieldMarker,13)
        %         keyboard;
        %         fieldHeader = strsplit(a);
        %         Ncomp = length(fieldHeader)-1;  % Count how many columns to expect from the header
        Ncomp = sum(ismember(a,'"'))/2; % Count how many columns to expect from the header
        getP = false;
        if Ncomp >=9, getP = true; end  % Assume column 9 is directivity..
        if fCount == 0
            % Initialise matrices
            NfInit = 1; % Use as initial maximum guess...
            [th,ph,Eth,Eph,Dth,Dph,Dtotclear] = deal(zeros(Nth*Nph,NfInit));
            fData = zeros(Nth*Nph,Ncomp,NfInit);
            Prad = zeros(1,NfInit);
        end
        fCount = fCount + 1;
        fDataForm = repmat('%f',1,Ncomp);
        fData(:,:,fCount) = fscanf(fid,fDataForm,[Ncomp,Nth*Nph])';
        th(:,fCount) = d2r(fData(:,1,fCount));
        ph(:,fCount) = d2r(fData(:,2,fCount));
        Eth(:,fCount) = fData(:,3,fCount) + 1i.*fData(:,4,fCount);
        Eph(:,fCount) = fData(:,5,fCount) + 1i.*fData(:,6,fCount);
        %         Dth(:,fCount) = fData(:,7,fCount);
        %         Dph(:,fCount) = fData(:,8,fCount);
        %         Dtot(:,fCount) = fData(:,9,fCount);
        % Get the power from the directivity/gain
        if getP
            U = 1/(2*eta0).*(abs(Eth(:,fCount)).^2 + abs(Eph(:,fCount)).^2);
            Prad(fCount) = median(4.*pi.*U./10.^(fData(:,9,fCount)./10));
        end
    end
end


%% Build the object

x = ph(:,1);
y = th(:,1);

E1 = Eth;
E2 = Eph;
E3 = zeros(size(E1)); %.ffe files will never have a radial field component...

coorSys = lower(coorSys); %coordinate system string fetched from header text of .ffe file
polType = 'linear'; %It seems that FEKO always outputs linear polarised fields (corresponding to th-ph coordinates)
gridType = 'PhTh'; %It seems that FEKO always outputs .ffe field values in theta-phi form, so this is hardcoded as such
%NB: Prad defined earlier
radEff = ones(size(freq)); %replace with manual radiation efficiency calculation
freqUnit = 'Hz'; %It seems that FEKO always outputs frequencies in Hz, so this is hardcoded as such

FF = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
% FF.polType = polType;
% FF.coorSys = coorSys;
% FF.gridType = gridType;
% FF.freqUnit = freqUnit;
FF = setEnames(FF);
FF = setXYnames(FF);
FF = setPhTh(FF);
FF = setFreq(FF);
FF = setBase(FF);


