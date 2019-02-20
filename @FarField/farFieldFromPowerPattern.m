function FF = farFieldFromPowerPattern(x,y,P,freq,Pdim,coorSys,polType,gridType,freqUnit,slant)

%doc goes here

%constants
load EMconstants
Nf = length(freq);

%Handle optional input arguments
optionNames = {'coorSys','polType','gridType','freqUnit','slant'};
defaultOptions = {'''spherical''','''linear''','''PhTh''','''Hz''','''0'''};
for ii = 1:length(optionNames) %run through the optional inputs, check if they have been provided, and if not set them to their defaults
    if (nargin < ii + 5)
        eval([optionNames{ii},' = ',defaultOptions{ii},';']);
    end
end

%NB: currently assuming that P is radiation intensity
E = (2*eta0).*sqrt(P); %Find E-field magnitude to distribute equally into both field components

%Check power pattern and generate E1 and E2 accordingly
switch Pdim
    case 1 % One-dimensional power pattern provided, generate axially symmetric pattern and E1, E2
        Nx = length(x)/length(P); %Find number of x-cuts, assuming P is given as a set of y-dependent x-cuts (e.g theta-dependent phi-cuts in spherical coords)
        if (Nx ~= floor(Nx)) error('Power pattern size does not fit an integer no. of times into x or y vectors'); end
        E1 = repmat(E,Nx,1)./2;
        E2 = E1;
    case 2 % Two-dimensional power pattern provided, generate corresponding E1, E2
        E1 = E./2; 
        E2 = E1; 
    otherwise
        error('Incompatible power pattern matrix provided - check dimensions of P')
end

%% Build the object

E3 = zeros(size(E1)); %.ffe files will never have a radial field component...

Prad = 4*pi.*ones(1,Nf); %!!! NB: this should be replaced with a 2D power integral of P
radEff = ones(size(freq));

FF = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);

FF = setEnames(FF);
FF = setXYnames(FF);
FF = setPhTh(FF);
FF = setFreq(FF);
FF = setBase(FF);
