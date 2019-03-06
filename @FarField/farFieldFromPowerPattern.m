function FF = farFieldFromPowerPattern(th,ph,P,freq,fieldPol)

%Name: farFieldFromPowerPattern.m
%Description:
%   Function to create a Farfield object from a power pattern, which must
%   be provided in the format that is outputted by powerPattern.m.
%Inputs:
% --th: column vector [Nang x 1] of th angles in rad
% --ph: column vector [Nang x 1] of ph angles in rad
% --P: column vector [Nang x 1] of power pattern values (see Pout output for powerPattern.m for an example of the format) 
% --freq: scalar frequency in Hz
% --fieldPol: (optional) string denoting how far field should be polarized- options are 'linearX', 'linearY', 'circularLH', 'circularRH' 
%Outputs:
% --FF: Farfield object containing parameters as determined from inputs


%constants
load EMconstants
Nf = length(freq);

%Handle optional input arguments
optionNames = {'fieldPol','coorSys','polType','gridType','freqUnit','slant'};
defaultOptions = {'''linearY''','''spherical''','''linear''','''PhTh''','''Hz''','''0'''};
for ii = 1:length(optionNames) %run through the optional inputs, check if they have been provided, and if not set them to their defaults
    if (nargin < ii + 4)
        eval([optionNames{ii},' = ',defaultOptions{ii},';']);
    end
end

%From power pattern and polarization parameters, generate E1 and E2 accordingly
switch fieldPol
    case 'linearX' % linearly polarised along X-axis 
        coorSys = 'Ludwig3';
        E1  = sqrt(P).*(2*eta0);
        E2  = zeros(size(P));
    case 'linearY' % linearly polarised along Y-axis 
        coorSys = 'Ludwig3';
        E1  = zeros(size(P));
        E2  = sqrt(P).*(2*eta0);
    case 'circularLH'  % Lefthand Circular polarization
        polType = 'circular';
        E1  = sqrt(P).*(2*eta0);
        E2  = zeros(size(P));
    case 'circularRH'  % Righthand Circular polarization
        polType = 'circular';
        E1  = zeros(size(P));
        E2  = sqrt(P).*(2*eta0);
    otherwise
        error('fieldPol input string unrecognised')
end
        


%% Build the object
E3 = zeros(size(E1)); %.ffe files will never have a radial field component...

Prad = 4*pi.*ones(1,Nf); %!!! NB: this should be replaced with a 2D power integral of P
radEff = ones(size(freq));

FF = FarField(ph,th,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);

FF = setEnames(FF);
FF = setXYnames(FF);
FF = setPhTh(FF);
FF = setFreq(FF);
FF = setBase(FF);