% function [] = plot(FF,freqIndex,plotType,varargin)
function [] = plot(FF,varargin)

% function plot(FF,name,value)
% Plots a representation of the farfield object in FF
% name, value are name value pairs, and can be the following:
%
% freqIndex is the index(es) of frequencies to be plotted (default 1)
%
% plotType can be:
%   ('3D') | '2D' | 'polar' | 'cartesian'
%
% output can be:
%   ('Directivity') | 'Gain' | 'E1' | 'E2' | 'AxialRatio' | 'AxialRatioInv'
%
% outputType can be:
%   ('mag') | 'phase' - only used for on E-field plots
%
% norm is a boolean (false) to normalize to maximum magnitude
%
% dynamicRange_dB is a (positive) dB value for the magnitude plot dynamic
% range (40)
%
% scaleMag can be:
%   ('dB') | 'lin' - only used for magnitude plots

% scalePhase can be:
%   ('deg') | 'rad' - only used for phase plots

narginchk(1,14);

%% Parsing through the inputs
parseobj = inputParser;
parseobj.FunctionName = 'plot';

typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'plot','FF',1);
addRequired(parseobj,'FF',typeValidationObj);

typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'vector','nonempty','integer'},'plot','freqIndex');
addParameter(parseobj,'freqIndex',1,typeValidationFreq);

typeValidationnorm = @(x) validateattributes(x,{'numeric'},{'binary','nonempty','numel',1},'plot','norm');
addParameter(parseobj,'norm',false,typeValidationnorm );

typeValidationDR = @(x) validateattributes(x,{'numeric'},{'real','positive','nonempty','numel',1},'plot','dynamicRange_dB');
addParameter(parseobj,'dynamicRange_dB',40,typeValidationDR );

expectedplotType = {'3D','2D','polar','cartesian'};
addParameter(parseobj,'plotType','3D', @(x) any(validatestring(x,expectedplotType)));

expectedoutput = {'Directivity','Gain','E1','E2','AxialRatio','AxialRatioInv'};
addParameter(parseobj,'output','Directivity', @(x) any(validatestring(x,expectedoutput)));

expectedoutputType = {'mag','phase'};
addParameter(parseobj,'outputType','mag', @(x) any(validatestring(x,expectedoutputType)));

expectedscaleMag = {'dB','lin'};
addParameter(parseobj,'scaleMag','dB', @(x) any(validatestring(x,expectedscaleMag)));

expectedscalePhase = {'deg','rad'};
addParameter(parseobj,'scalePhase','deg', @(x) any(validatestring(x,expectedscalePhase)));

parse(parseobj, FF, varargin{:});

%% Extract plot output type
freqIndex = parseobj.Results.freqIndex;
output = parseobj.Results.output;
outputType = parseobj.Results.outputType;
norm = parseobj.Results.norm;
dynamicRange_dB = parseobj.Results.dynamicRange_dB;
plotType = parseobj.Results.plotType;
scaleMag = parseobj.Results.scaleMag;
scalePhase = parseobj.Results.scalePhase;

switch output
    case {'Directivity','Gain','AxialRatio','AxialRatioInv'} 
        dr = lin10(-dynamicRange_dB);
        if strcmp(output,'Directivity')
            Zmat = getDirectivity(FF);
            compName = 'Directivity';
        elseif strcmp(output,'Gain')
            Zmat = getGain(FF);
            compName = 'Gain (IEEE)';
        elseif strcmp(output,'AxialRatio')
            [Zmat,~] = getAxialRatio(FF);
            compName = 'Axial Ratio';
        elseif strcmp(output,'AxialRatioInv')
            [~,Zmat] = getAxialRatio(FF);
            compName = 'Inverse Axial Ratio';
        end
        if norm, Zmat = bsxfun(@times,Zmat,1./(max(Zmat))); end
        if strcmp(scaleMag,'dB')
            Zmat = dB10(Zmat); unit = 'dB';
            if norm, unit = 'dBi'; end
        else
            unit = 'pu';
        end
    case {'E1','E2'}
        dr = lin20(-dynamicRange_dB);
        if strcmp(output,'E1')
            [Zmat,~,~] = getEfield(FF);
            compName = FF.E1name;
        elseif strcmp(output,'E2')
            [~,Zmat,~] = getEfield(FF);
            compName = FF.E2name;
        end
        if strcmp(outputType,'mag')
            Zmat = abs(Zmat);
            if norm, Zmat = bsxfun(@times,Zmat,1./(max(Zmat))); end
            if strcmp(scaleMag,'dB')
                Zmat = dB20(Zmat); unit = 'dBV/m';
                if norm, unit = 'dB'; end
            else
                if norm, unit = ''; else unit = 'V/m'; end  
            end
        elseif strcmp(outputType,'phase')
            Zmat = angle(Zmat);
            if norm, Zmat = bsxfun(@plus,Zmat,1./(min(Zmat))); end
            unit = 'rad';
            if strcmp(scalePhase,'deg'), Zmat = rad2deg(Zmat); unit = 'deg'; end
        end
end

Z = Zmat(:,freqIndex);
freqVect = FF.freq(freqIndex);

%% Unpack angles
th_vect = unique(FF.th);
ph_vect = unique(FF.ph);
[PH,TH] = meshgrid(ph_vect,th_vect);


%% Plot for each frequency
for ff = 1:numel(freqVect)
    Zplot = Z(:,ff);
    % Fix Dynamic Range if required
    if strcmp(outputType,'mag')
        if strcmp(scaleMag,'dB')
            Zplot(Zplot < (max(Zplot) - dynamicRange_dB)) = max(Zplot) - dynamicRange_dB;
        else
            Zplot(Zplot < max(Zplot)*dr) = max(Zplot)*dr;
        end
    end
    
    figure
    switch plotType
        case '3D'
            % Use the MATLAB antennas toolbox plotting function
            patternCustom(Zplot,rad2deg(FF.th),rad2deg(FF.ph));
        case '2D'
            surf(rad2deg(TH),rad2deg(PH),reshape(Zplot,FF.Nth,FF.Nph))
            xlabel('\theta (deg)')
            ylabel('\phi (deg)')
    end
    
%     if ~magFlag
%         mp_string = 'phase';
%     end
    title([FF.polBase, ', ',FF.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,')'])

end
% keyboard;








end