% function [] = plot(FF,freqIndex,plotType,varargin)
function [] = plot(FF,varargin)

% function plot(FF,name,value)
% Plots a representation of the farfield object in FF
% name, value are name value pairs, and can be the following:
%
% freqIndex is the index of the frequency to be plotted (default 1)
%
% plotType can be:
%   ('3D') | '2D' | 'polar' | 'cartesian'
%
% output can be:
%   ('Directivity') | 'Gain' | 'E1' | 'E2' | 'AxialRatio' | 'AxialRatioInv'
%   'CO/XP' | 'XP/CO'
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
%
% freqUnit can be:
%   ('GHz') | 'Hz' | 'kHz' | 'MHz' | 'THz'
%
% projection can be:
%   ('std') | 'TrueView' | 'DirCosine' | 'ArcSin'
%    See Masters and Gregson paper for details.  'std' just uses the
%    standard angle definitions of the polBase
%
% cutConstant can be (used only for polar and cartesian plots):
%   ('x') | 'y'  (x = [ph|az|ep|u|Xg|asin(u)]; y = [th|el|al|Yg|asin(v)]
%
% cutValue can be any value in the available angle range.  If omitted, the
% principle plane cuts will be attempted...
%
% step is the plot step size.  Can be empty - then the available data will
% be used.  If not, a griddata interpolant will be made.
%
% plotProperties can be a variety of name, value pairs including:
%   LineWidth, LineStyle (like '-k.')

narginchk(1,20);

%% Parsing through the inputs
parseobj = inputParser;
parseobj.FunctionName = 'plot';

typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'plot','FF',1);
addRequired(parseobj,'FF',typeValidationObj);

typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'plot','freqIndex');
addParameter(parseobj,'freqIndex',1,typeValidationFreq);

typeValidationnorm = @(x) validateattributes(x,{'numeric'},{'binary','nonempty','numel',1},'plot','norm');
addParameter(parseobj,'norm',false,typeValidationnorm );

typeValidationDR = @(x) validateattributes(x,{'numeric'},{'real','positive','nonempty','numel',1},'plot','dynamicRange_dB');
addParameter(parseobj,'dynamicRange_dB',40,typeValidationDR );

expectedplotType = {'3D','2D','polar','cartesian'};
addParameter(parseobj,'plotType','3D', @(x) any(validatestring(x,expectedplotType)));

expectedoutput = {'Directivity','Gain','E1','E2','AxialRatio','AxialRatioInv','CO/XP','XP/CO'};
addParameter(parseobj,'output','Directivity', @(x) any(validatestring(x,expectedoutput)));

expectedoutputType = {'mag','phase'};
addParameter(parseobj,'outputType','mag', @(x) any(validatestring(x,expectedoutputType)));

expectedscaleMag = {'dB','lin'};
addParameter(parseobj,'scaleMag','dB', @(x) any(validatestring(x,expectedscaleMag)));

expectedscalePhase = {'deg','rad'};
addParameter(parseobj,'scalePhase','deg', @(x) any(validatestring(x,expectedscalePhase)));

expectedfreqUnit = {'Hz','kHz','MHz','GHz','THz'};
addParameter(parseobj,'freqUnit','GHz', @(x) any(validatestring(x,expectedfreqUnit)));

expectedprojection = {'std','TrueView','DirCosine','ArcSin'};
addParameter(parseobj,'projection','std', @(x) any(validatestring(x,expectedprojection)));

expectedcutConstant = {'x','y'};
addParameter(parseobj,'cutConstant','x', @(x) any(validatestring(x,expectedcutConstant)));

typeValidationcutValue = @(x) validateattributes(x,{'numeric'},{'real'},'plot','cutValue');
addParameter(parseobj,'cutValue',[],typeValidationcutValue);

typeValidationstep = @(x) validateattributes(x,{'numeric'},{'real'},'plot','step');
addParameter(parseobj,'step',[],typeValidationstep);

typeValidationLineWidth = @(x) validateattributes(x,{'numeric'},{'real'},'plot','LineWidth');
addParameter(parseobj,'LineWidth',1,typeValidationLineWidth);

addParameter(parseobj,'LineStyle',[]);

parse(parseobj, FF, varargin{:});

freqIndex = parseobj.Results.freqIndex;
output = parseobj.Results.output;
outputType = parseobj.Results.outputType;
norm = parseobj.Results.norm;
dynamicRange_dB = parseobj.Results.dynamicRange_dB;
plotType = parseobj.Results.plotType;
scaleMag = parseobj.Results.scaleMag;
scalePhase = parseobj.Results.scalePhase;
freqUnit = parseobj.Results.freqUnit;
projection = parseobj.Results.projection;
cutConstant = parseobj.Results.cutConstant;
cutValue = parseobj.Results.cutValue;
step = parseobj.Results.step;
LineWidth = parseobj.Results.LineWidth;
LineStyle = parseobj.Results.LineStyle;


%% Extract plot output type

switch output
    case {'Directivity','Gain','AxialRatio','AxialRatioInv','CO/XP','XP/CO'}
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
        elseif strcmp(output,'CO/XP') || strcmp(output,'XP/CO')
            switch output
                case 'CO/XP'
                    num = abs(FF.E2);
                    den = abs(FF.E1);
                    compName = 'Co-pol/Cross-pol';
                case 'XP/CO'
                    num = abs(FF.E1);
                    den = abs(FF.E2);
                    compName = 'Cross-pol/Co-pol';
            end
            % Sort out dynamic range here to avoid divide by zero
            maxVal = max(max(num));
            minVal = lin20(-dynamicRange_dB).*maxVal;
            den(den < minVal) = minVal;
            Zmat = (num./den).^2;
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
% freqVect = FF.freq(freqIndex);

%% Extract some angles
switch projection
    case 'std'
        switch FF.polBase
            case {'spherical','Ludwig1','Ludwig3'}
                x = rad2deg(FF.ph);
                y = rad2deg(FF.th);
                xname = '\phi (deg)';
                yname = '\theta (deg)';
            case 'Ludwig2AE'
                [y,x] = getElAz(FF);
                x = rad2deg(x);
                y = rad2deg(y);
                xname = 'Azimuth (deg)';
                yname = 'Elevation (deg)';
            case 'Ludwig2EA'
                [y,x] = getElAz(FF);
                x = rad2deg(x);
                y = rad2deg(y);
                xname = 'Epsilon (deg)';
                yname = 'Alpha (deg)';
        end
    case 'TrueView'
        [x,y] = getXgYg(FF);
        x = rad2deg(x);
        y = rad2deg(y);
        xname = 'X_g = \theta cos(\phi) (deg)';
        yname = 'Y_g = \theta sin(\phi) (deg)';
    case 'DirCosine'
        [x,y,~] = getUVW(FF);
        xname = 'u = sin(\theta)cos(\phi)';
        yname = 'v = sin(\theta)sin(\phi)';
    case 'ArcSin'
        [u,v,~] = getUVW(FF);
        x = rad2deg(asin(u));
        y = rad2deg(asin(v));
        xname = 'asin(u) (deg)';
        yname = 'asin(v) (deg)';
end

%% Make the plots
switch freqUnit
    case 'Hz'
        freqMult = 1;
    case 'kHz'
        freqMult = 1e-3;
    case 'MHz'
        freqMult = 1e-6;
    case 'GHz'
        freqMult = 1e-9;
    case 'THz'
        freqMult = 1e-12;
end
freqPlot = FF.freq(freqIndex)*freqMult;
Zplot = Z(:,freqIndex);
% Fix Dynamic Range if required
if strcmp(outputType,'mag')
    if strcmp(scaleMag,'dB')
        Zplot(Zplot < (max(Zplot) - dynamicRange_dB)) = max(Zplot) - dynamicRange_dB;
    else
        Zplot(Zplot < max(Zplot)*dr) = max(Zplot)*dr;
    end
end

xMat = reshape(x,FF.Nth,FF.Nph);
yMat = reshape(y,FF.Nth,FF.Nph);
if isempty(step)
    X = xMat;
    Y = yMat;
elseif ~strcmp(plotType,'3D') && ~isempty(cutValue)
    xvect = min(x):step:max(x);
    yvect = min(y):step:max(y);
    [X,Y] = meshgrid(xvect,yvect);
    Zplot = griddata(xMat,yMat,reshape(Zplot,FF.Nth,FF.Nph),X,Y);
    Zplot = Zplot(:);
    x = X(:);
    y = Y(:);
end

%     figure
switch plotType
    case '3D'
        % Use the MATLAB antennas toolbox plotting function
        patternCustom(Zplot,rad2deg(FF.th),rad2deg(FF.ph));
        title([FF.polBase, ', ',FF.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,'); Freq = ',num2str(freqPlot),' ', freqUnit])
    case '2D'
%         surf(X,Y,reshape(Zplot,FF.Nth,FF.Nph),'EdgeColor','Interp','FaceColor','Interp')
        surf(X,Y,reshape(Zplot,size(X)),'EdgeColor','Interp','FaceColor','Interp')
        xlabel(xname)
        ylabel(yname)
        view([0,90])
        axis equal
        xlim([min(x),max(x)])
        ylim([min(y), max(y)])
        colorbar
        title([FF.polBase, ', ',FF.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,'); Freq = ',num2str(freqPlot),' ', freqUnit])
    case 'cartesian'
        lw = LineWidth;   % Can maybe chance to be accessed form outside later...
        if isempty(cutValue)
            % Find the principle cuts
            iph0 = find(abs(FF.ph - 0) < eps);
            iph45 = find(abs(FF.ph - deg2rad(45)) < eps);
            iph90 = find(abs(FF.ph - deg2rad(90)) < eps);
            iph135 = find(abs(FF.ph - deg2rad(135)) < eps);
            iph180 = find(abs(FF.ph - deg2rad(180)) < eps);
            iph225 = find(abs(FF.ph - deg2rad(225)) < eps);
            iph270 = find(abs(FF.ph - deg2rad(270)) < eps);
            iph315 = find(abs(FF.ph - deg2rad(315)) < eps);
            
            plot(rad2deg(FF.th(iph0)),Zplot(iph0),[LineStyle,'k'],'lineWidth',lw), grid on, hold on
            plot(rad2deg(FF.th(iph45)),Zplot(iph45),[LineStyle,'b'],'lineWidth',lw)
            plot(rad2deg(FF.th(iph90)),Zplot(iph90),[LineStyle,'r'],'lineWidth',lw)
            plot(rad2deg(FF.th(iph135)),Zplot(iph135),[LineStyle,'g'],'lineWidth',lw)
            plot(-rad2deg(FF.th(iph180)),Zplot(iph180),[LineStyle,'k'],'lineWidth',lw)
            plot(-rad2deg(FF.th(iph225)),Zplot(iph225),[LineStyle,'b'],'lineWidth',lw)
            plot(-rad2deg(FF.th(iph270)),Zplot(iph270),[LineStyle,'r'],'lineWidth',lw)
            plot(-rad2deg(FF.th(iph315)),Zplot(iph315),[LineStyle,'g'],'lineWidth',lw)
            lg = legend('\phi = 0^\circ','\phi = 45^\circ','\phi = 90^\circ','\phi = 135^\circ');
            xlabel('\theta (deg)')
            ylabel([outputType,'(', compName, ') (',unit,')'])
            %                 lg.String = lg.String(1:4); % Just keep the first four legend entries
        else
            if cutConstant == 'x'
                iCut = find(abs(x - cutValue) < eps);
                plot(y(iCut),Zplot(iCut),LineStyle,'lineWidth',lw), grid on, hold on
                xlabel(yname)
            elseif cutConstant == 'y'
                iCut = find(abs(y - cutValue) < eps);
                plot(x(iCut),Zplot(iCut),LineStyle,'lineWidth',lw), grid on, hold on
                xlabel(xname)
            end
            ylabel([outputType,'(', compName, ') (',unit,')'])
            % ToDo
        end
        title([FF.polBase, ', ',FF.polType, ' polarisation; Freq = ',num2str(freqPlot),' ', freqUnit])

end


% end

end