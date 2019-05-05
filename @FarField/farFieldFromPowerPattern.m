function FF = farFieldFromPowerPattern(x,y,P,freq,fieldPol,coorSys,gridType,freqUnit,slant)

%Name: farFieldFromPowerPattern.m
%Description:
%   Function to create a Farfield object from a power pattern, which must
%   be provided in the format that is outputted by powerPattern.m.
%Inputs:
% --x: column vector [Nang x 1] of x-axis angular coordinates in rad
% --y: column vector [Nang x 1] of y-axis in rad
% --P: column vector [Nang x Nf] of power pattern values (see Pout output for powerPattern.m for an example of the format) 
% --freq: scalar frequency in Hz
% --fieldPol: (optional) string denoting how far field should be polarized- options are 'linearX', 'linearY', 'circularLH', 'circularRH' 
%Outputs:
% --FF: Farfield object containing parameters as determined from inputs


%constants
load EMconstants
Nf = length(freq);

%Handle optional input arguments
optionNames = {'fieldPol','coorSys','gridType','polType','freqUnit','slant'};
defaultOptions = {'''linearY''','''spherical''','''PhTh''','''linear''','''Hz''','''0'''};
for ii = 1:length(optionNames) %run through the optional inputs, check if they have been provided, and if not set them to their defaults
    if (nargin < ii + 4) || isempty(eval(optionNames{ii}))
        eval([optionNames{ii},' = ',defaultOptions{ii},';']);
    end
end
if isscalar(P)
    P = P./(4.*pi).*ones(size(x));
end
%From power pattern and polarization parameters, generate E1 and E2 accordingly
switch fieldPol
    case 'linearX' % linearly polarised along X-axis 
        coorSys = 'Ludwig3';
        E1  = sqrt(P.*2*eta0);
        E2  = zeros(size(P));
    case 'linearY' % linearly polarised along Y-axis 
        coorSys = 'Ludwig3';
        E1  = zeros(size(P));
        E2  = sqrt(P.*2*eta0);
    case 'circularLH'  % Lefthand Circular polarization
        polType = 'circular';
        E1  = sqrt(P.*2*eta0);
        E2  = zeros(size(P));
    case 'circularRH'  % Righthand Circular polarization
        polType = 'circular';
        E1  = zeros(size(P));
        E2  = sqrt(P.*2*eta0);
    otherwise
        error('fieldPol input string unrecognised')
end
        
%% Build the object
E3 = zeros(size(E1)); %farfields will never have a radial field component...

Prad = ones(1,Nf); % Dummy for the constructor
radEff = ones(size(freq));

FF = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);

FF.Prad = FF.pradInt;
FF.Directivity_dBi = dB10(max(FF.getDirectivity()));
FF.Gain_dB = dB10(max(FF.getGain()));
FF = setEnames(FF);
FF = setXYnames(FF);
FF = setPhTh(FF);
FF = setFreq(FF);
FF = setBase(FF);
