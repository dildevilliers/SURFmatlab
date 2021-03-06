function FF = farFieldFromPowerPattern(ph,th,P,freq,fieldPol,freqUnit)

%Name: farFieldFromPowerPattern.m
%Description:
%   Function to create a Farfield object from a power pattern, which must
%   be provided in the format that is outputted by powerPattern.m.
%Inputs:
% --ph: column vector [Nang x 1] of ph angles in rad
% --th: column vector [Nang x 1] of th angles in rad
% --P: column vector [Nang x 1] of power pattern values (see Pout output for powerPattern.m for an example of the format) 
% --freq: scalar frequency in Hz
% --fieldPol: (optional) string denoting how far field should be polarized- options are 'linearX', 'linearY', 'circularLH', 'circularRH' 
% --freqUnit: (optional) string denoting the frequency unit {('Hz'),'kHz','MHz','GHz','THz'}
%Outputs:
% --FF: Farfield object containing parameters as determined from inputs


%constants
load EMconstants
Nf = length(freq);

%Handle optional input arguments
optionNames = {'fieldPol','freqUnit'};
defaultOptions = {'''linearY''','''Hz'''};
for ii = 1:length(optionNames) %run through the optional inputs, check if they have been provided, and if not set them to their defaults
    if (nargin < ii + 4)
        eval([optionNames{ii},' = ',defaultOptions{ii},';']);
    end
end
if isscalar(P)
    P = P./(4.*pi).*ones(size(ph));
end
%From power pattern and polarization parameters, generate E1 and E2 accordingly
coorType = 'Ludwig3';
gridType = 'PhTh';
switch fieldPol
    case 'linearX' % linearly polarised along X-axis 
        polType = 'linear';
        E1  = sqrt(P.*2*eta0);
        E2  = zeros(size(P));
    case 'linearY' % linearly polarised along Y-axis 
        polType = 'linear';
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
E3 = [];
Prad = ones(1,Nf); % Dummy for the constructor
radEff = ones(size(freq));

FF = FarField(ph,th,E1,E2,E3,freq,Prad,radEff,coorType,polType,gridType,freqUnit);
FF.Prad = FF.pradInt;

end