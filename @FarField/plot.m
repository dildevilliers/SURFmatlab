function [] = plot(obj,varargin)

% function plot(obj,name,value)
% Plots a representation of the farfield object in obj
% name, value are name value pairs, and can be the following:
%
% freqIndex is the index of the frequency to be plotted (default 1)
%
% plotType can be:
%   ('3D') | '2D' | 'polar' | 'cartesian'
%
% output can be:
%   ('Directivity') | 'Gain' | 'E1' | 'E2' | 'AxialRatio' | 'AxialRatioInv'
%   'CO_XP' | 'XP_CO' | 'W' | 'U'
%
% outputType can be:
%   ('mag') | 'phase' | 'real' | 'imag' last 3 only used for on E-field plots
%
% norm is a boolean (false) to normalize to maximum magnitude
%
% dynamicRange_dB is a (positive) dB value for the magnitude plot dynamic
% range (40)
%
% scaleMag can be:
%   ('dB') | 'lin' - only used for magnitude plots
%
% scalePhase can be:
%   ('deg') | 'rad' - only used for phase plots
%
% freqUnit can be:
%   ('GHz') | 'Hz' | 'kHz' | 'MHz' | 'THz'
%
% cutConstant can be (used only for polar and cartesian plots):
%   ('x') | 'y'  (x = [ph|az|ep|u|Xg|asin(u)]; y = [th|el|al|Yg|asin(v)]
%
% cutValue can be any value in the available angle range (in rad).  
%
% step is the plot step size.  Can be empty - then the available data will
% be used and no surface will be plotted.  If not, a griddata interpolant will be made.
%
% plotProperties can be a variety of name, value pairs including:
%   LineWidth, LineStyle, Color (like '-.')
%
% showGrid is a boolean (false) to show the 2D grid where the data is
% calculated before interpolation.
%
% hemisphere is used in gridTypes DirCos and ArcSin and can be:
%   ('top') | 'bot'



narginchk(1,40);

%% Parsing through the inputs
parseobj = inputParser;
parseobj.FunctionName = 'plot';

typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'plot','obj',1);
addRequired(parseobj,'obj',typeValidationObj);

typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'plot','freqIndex');
addParameter(parseobj,'freqIndex',1,typeValidationFreq);

typeValidationnorm = @(x) validateattributes(x,{'logical','numeric'},{'binary','nonempty','numel',1},'plot','norm');
addParameter(parseobj,'norm',false,typeValidationnorm );

typeValidationDR = @(x) validateattributes(x,{'numeric'},{'real','positive','nonempty','numel',1},'plot','dynamicRange_dB');
addParameter(parseobj,'dynamicRange_dB',40,typeValidationDR );

expectedplotType = {'3D','2D','polar','cartesian'};
addParameter(parseobj,'plotType','3D', @(x) any(validatestring(x,expectedplotType)));

expectedoutput = {'Directivity','Gain','E1','E2','E3','AxialRatio','AxialRatioInv','CO_XP','XP_CO','W','U'};
addParameter(parseobj,'output','Directivity', @(x) any(validatestring(x,expectedoutput)));

expectedoutputType = {'mag','phase','real','imag'};
addParameter(parseobj,'outputType','mag', @(x) any(validatestring(x,expectedoutputType)));

expectedscaleMag = {'dB','lin'};
addParameter(parseobj,'scaleMag','dB', @(x) any(validatestring(x,expectedscaleMag)));

expectedscalePhase = {'deg','rad'};
addParameter(parseobj,'scalePhase','deg', @(x) any(validatestring(x,expectedscalePhase)));

expectedfreqUnit = {'Hz','kHz','MHz','GHz','THz'};
addParameter(parseobj,'freqUnit','GHz', @(x) any(validatestring(x,expectedfreqUnit)));

expectedcutConstant = {'x','y'};
addParameter(parseobj,'cutConstant','x', @(x) any(validatestring(x,expectedcutConstant)));

typeValidationcutValue = @(x) validateattributes(x,{'numeric'},{'real'},'plot','cutValue');
addParameter(parseobj,'cutValue',[],typeValidationcutValue);

typeValidationstep = @(x) validateattributes(x,{'numeric'},{'real'},'plot','step');
addParameter(parseobj,'step',[],typeValidationstep);     % In degrees

typeValidationLineWidth = @(x) validateattributes(x,{'numeric'},{'real'},'plot','LineWidth');
addParameter(parseobj,'LineWidth',1,typeValidationLineWidth);

typeValidationLineStyle = @(x) validateattributes(x,{'char'},{'nonempty'},'plot','LineStyle');
addParameter(parseobj,'LineStyle','-',typeValidationLineStyle);     

typeValidationColor = @(x) validateattributes(x,{'numeric','char'},{'nonempty'},'plot','Color');
addParameter(parseobj,'Color','k',typeValidationColor);     

typeValidationshowGrid = @(x) validateattributes(x,{'logical','numeric'},{'binary','nonempty','numel',1},'plot','showGrid');
addParameter(parseobj,'showGrid',false,typeValidationshowGrid );

expectedhemisphere = {'top','bot'};
addParameter(parseobj,'hemisphere','top', @(x) any(validatestring(x,expectedhemisphere)));

parse(parseobj, obj, varargin{:});

freqIndex = parseobj.Results.freqIndex;
output = parseobj.Results.output;
outputType = parseobj.Results.outputType;
norm = parseobj.Results.norm;
dynamicRange_dB = parseobj.Results.dynamicRange_dB;
plotType = parseobj.Results.plotType;
scaleMag = parseobj.Results.scaleMag;
scalePhase = parseobj.Results.scalePhase;
freqUnit = parseobj.Results.freqUnit;
cutConstant = parseobj.Results.cutConstant;
cutValue = parseobj.Results.cutValue;
step = parseobj.Results.step;
LineWidth = parseobj.Results.LineWidth;
LineStyle = parseobj.Results.LineStyle;
Color = parseobj.Results.Color;
showGrid = parseobj.Results.showGrid;
hemisphere = parseobj.Results.hemisphere;

%% Sort out the plot grid and names

% Get valid positions for the plot
if strcmp(obj.gridType,'DirCos') || strcmp(obj.gridType,'ArcSin')
    % Try to get the direction cosines from the base grid definition - if the
    % base definition is not a direction cosine type it can contain
    % information over the full sphere. 
    objBase = obj.grid2Base;
    grid2DirCoshandle = str2func([objBase.gridType,'2DirCos']);
    [~,~,w] = grid2DirCoshandle(objBase.x,objBase.y);
    if strcmp(hemisphere,'top')
        valAng = w >= 0;
    elseif strcmp(hemisphere,'bot')
        valAng = w <= 0;
    end
else
    valAng = ones(obj.Nang,1);
end

% Get the original grid and output
X = reshape(obj.x,obj.Ny,obj.Nx);
Y = reshape(obj.y,obj.Ny,obj.Nx);
if strcmp(output,'E1')
    [Z,~,~] = getEfield(obj);
elseif strcmp(output,'E2')
    [~,Z,~] = getEfield(obj);
elseif strcmp(output,'E3')
    [~,~,Z] = getEfield(obj);
else
    outputHandle = str2func(['get',output]);
    Z = outputHandle(obj);
end
Z = Z(:,freqIndex);

if isempty(step)
    if strcmp(plotType,'cartesian') || strcmp(plotType,'cartesian')
        error('Cant do 1D plot with empty step size - this is reserved for making 2D grids etc.')
    else
        Xi = X;
        Yi = Y;
        xi = obj.x;
        yi = obj.y;
        NxPlot = obj.Nx;
        NyPlot = obj.Ny;
        Zi = Z;
    end
else
    % Get the interpolated plot points from the step information
    if strcmp(obj.gridType,'DirCos') || strcmp(obj.gridType,'ArcSin')
        step = sind(step);
    else
        step = deg2rad(step);
    end
    xivect = min(obj.x):step:max(obj.x);
    yivect = min(obj.y):step:max(obj.y);
    [Xi,Yi] = meshgrid(xivect,yivect);
    NxPlot = numel(xivect);
    NyPlot = numel(yivect);
    switch plotType
        case {'3D','2D'}
            xi = Xi(:);
            yi = Yi(:);
        case {'cartesian','polar'}
            switch cutConstant
                case 'x'
                    yi = yivect;
                    xi = ones(size(yi)).*cutValue;
                case 'y'
                    xi = xivect;
                    yi = ones(size(xi)).*cutValue;
            end
    end
    [Zi] = interpolateGrid(obj,output,xi,yi,freqIndex,hemisphere);
end

% Assign axis names
switch obj.gridType
    case 'DirCos'
        xiplot = xi;
        yiplot = yi;
        xname = [obj.xname, ' = sin(\theta)cos(\phi)'];
        yname = [obj.yname, ' = sin(\theta)sin(\phi)'];
        axisUnit = '';
    otherwise
        X = rad2deg(X);
        Y = rad2deg(Y);
        Xi = rad2deg(Xi);
        Yi = rad2deg(Yi);
        xiplot = rad2deg(xi);
        yiplot = rad2deg(yi);
        axisUnit = '(deg)';
        xname = [obj.xname, ' ' ,axisUnit];
        yname = [obj.yname, ' ' ,axisUnit];
end

%% Condition outputs
% Phase results
if (strcmp(output,'E1') || strcmp(output,'E2') || strcmp(output,'E3')) && strcmp(outputType,'phase')
            Zplot = angle(Z);
            Zplot(~valAng) = NaN;
            Ziplot = angle(Zi);
            if norm
                Zplot = Zplot - max(Zplot);
                Ziplot = Ziplot - max(Ziplot);
            end
            unit = 'V/m (rad)';
            if strcmp(scalePhase,'deg')
                Zplot = rad2deg(Zplot);
                Ziplot = rad2deg(Ziplot);
                unit = 'V/m (deg)';
            end
            compName = strrep(obj.([output,'name']),'_','\');
else    % Magnitude/Component results
    switch outputType
        case 'mag'
            Zplot = abs(Z);
            Ziplot = abs(Zi);
        case 'real'
            Zplot = real(Z);
            Ziplot = real(Zi);
        case 'imag'
            Zplot = imag(Z);
            Ziplot = imag(Zi);
    end
    Zplot(~valAng) = NaN;
    switch output
        case {'Directivity','Gain','AxialRatio','AxialRatioInv','CO_XP','XP_CO','W','U'}
            dBscale = 10;
            unit = '';
            compName = strrep(output,'_','/');
        case {'E1','E2','E3'}
            dBscale = 20;
            unit = 'V/m ';
            compName = strrep(obj.([output,'name']),'_','\');
    end
    if norm
        Zplot = Zplot./max(Zplot);
        Ziplot = Ziplot./max(Ziplot);
    end
    if strcmp(scaleMag,'dB')
        dBHandle = str2func(['dB',num2str(dBscale)]);
        Zplot = dBHandle(Zplot);
        Ziplot = dBHandle(Ziplot);
        unit = [unit, 'dB'];
    end
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
freqPlot = obj.freqHz(freqIndex)*freqMult;

switch plotType
    case '3D'
        % ToDo:
        % Doesn't work with AzEl/ElAz grids +-90 y-axis breaks the plot
        
        % Use the MATLAB antennas toolbox plotting function
        if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
            % Handle dynamic range here: ToDo
            if strcmp(outputType,'mag')
                maxVal = max(Zplot(:));
                switch scaleMag
                    case 'dB'
                        if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                            maxVal = (1-norm).*dynamicRange_dB;
                            minVal = norm.*dynamicRange_dB;
                        else
                            minVal = maxVal-dynamicRange_dB;
                        end
                        Ziplot(Ziplot<minVal) = minVal;
                        Ziplot(Ziplot>maxVal) = maxVal;
                    case 'lin'
                        linHandle = str2func(['lin',num2str(dBscale)]);
                        if ~norm
                            dr = linHandle(dynamicRange_dB);
                            if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                                caxis([0,dr]);
                            else
                                caxis([maxVal/dr,maxVal]);
                            end
                        end
                end
            end
            iVal = ~isnan(Ziplot);
            patternCustom(Ziplot(iVal),Yi(iVal),Xi(iVal));
            title([obj.coorSys, ', ',obj.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,'); Freq = ',num2str(freqPlot),' ', freqUnit])
        else
            error(['gridType must be PhTh for 3D plots: found gridType = ', obj.gridType])
        end
        
    case '2D'
        Ziplot = reshape(Ziplot,NyPlot,NxPlot);
        if ~isempty(step)
            surf(Xi,Yi,Ziplot,'EdgeColor','Interp','FaceColor','Interp')
            colorbar
        end
        if showGrid || isempty(step)
            hold on
            Zplot(valAng == 0) = NaN;
            plot3(X(:),Y(:),Zplot(:),'k.')
            hold off
        end
        xlabel(xname)
        ylabel(yname)
        view([0,90])
        axis equal
        xlim([min(xiplot),max(xiplot)])
        ylim([min(yiplot),max(yiplot)])
        % Handle dynamic range here
        if strcmp(outputType,'mag')
            maxVal = max(Zplot(:));
            switch scaleMag
                case 'dB'
                    if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                        rangeZ = [0,dynamicRange_dB] - norm.*dynamicRange_dB;
                        caxis(rangeZ);
%                         zlim(rangeZ);
                    else
                        rangeZ = [maxVal-dynamicRange_dB,maxVal];
                        caxis(rangeZ);
%                         zlim(rangeZ);
                    end
                case 'lin'
                    linHandle = str2func(['lin',num2str(dBscale)]);
                    if ~norm
                        dr = linHandle(dynamicRange_dB);
                        if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                            rangeZ = [0,dr];
                            caxis(rangeZ);
%                             zlim(rangeZ);
                        else
                            rangeZ = [maxVal/dr,maxVal];
                            caxis(rangeZ);
%                             zlim([rangeZ]);
                        end
                    end
            end
        end
        title([obj.coorSys, ', ',obj.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,'); Freq = ',num2str(freqPlot),' ', freqUnit])
    case {'cartesian','polar'}
        % Initial bookkeeping to seperate the two options
        if strcmp(plotType,'cartesian')
            plotHandle = str2func(['plot']);
            limitHandle = str2func(['ylim']);
            xscale = 1;
        elseif strcmp(plotType,'polar')
            assert(~strcmp(obj.gridType,'DirCos'),'Polar plots not supported for DirCos gridType');
            plotHandle = str2func(['polarplot']);
            limitHandle = str2func(['rlim']);
            xscale = pi/180;
        end
        switch cutConstant
            case 'x'
                plotHandle(yiplot.*xscale,Ziplot,'LineStyle',LineStyle,'LineWidth',LineWidth,'Color',Color), grid on
                xlab = yname;
                cutName = obj.xname;
            case 'y'
                plotHandle(xiplot.*xscale,Ziplot,'LineStyle',LineStyle,'LineWidth',LineWidth,'Color',Color), grid on
                xlab = xname;
                cutName = obj.yname;
        end
        % Handle dynamic range here
        if strcmp(outputType,'mag')
            maxVal = max(Zplot(:));
            switch scaleMag
                case 'dB'
                    if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                        limitHandle([0,dynamicRange_dB] - norm.*dynamicRange_dB);
                    else
                        limitHandle([maxVal-dynamicRange_dB,maxVal]);
                    end
                case 'lin'
                    linHandle = str2func(['lin',num2str(dBscale)]);
                    if ~norm
                        dr = linHandle(dynamicRange_dB);
                        if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                            limitHandle([0,dr]);
                        else
                            limitHandle([maxVal/dr,maxVal]);
                        end
                    end
            end
        end
        if ~strcmp(obj.gridType,'DirCos')
            cutValue = rad2deg(cutValue);
        end
        titText = [obj.coorSys, ', ',obj.polType, ' polarisation; Freq = ',num2str(freqPlot),' ', freqUnit,'; ',cutName, ' = ',num2str(cutValue), ' ',axisUnit];

        % Final bookkeeping to seperate the two options
        hold on
        ax = gca;
        ylab = [outputType,'(', compName, ') (',unit,')'];
        if strcmp(plotType,'cartesian')
            xlabel(xlab)
            ylabel(ylab)
            title(titText)
        elseif strcmp(plotType,'polar')
            title([ylab,'; ',titText])
            ax.ThetaZeroLocation = 'top';
            ax.ThetaDir = 'clockwise';
        end
        
%     case 'cartesian'
%         lw = LineWidth;   
%         if isempty(cutValue)
%             % Find the principle cuts
%             iph0 = find(abs(obj.ph - 0) < eps);
%             iph45 = find(abs(obj.ph - deg2rad(45)) < eps);
%             iph90 = find(abs(obj.ph - deg2rad(90)) < eps);
%             iph135 = find(abs(obj.ph - deg2rad(135)) < eps);
%             iph180 = find(abs(obj.ph - deg2rad(180)) < eps);
%             iph225 = find(abs(obj.ph - deg2rad(225)) < eps);
%             iph270 = find(abs(obj.ph - deg2rad(270)) < eps);
%             iph315 = find(abs(obj.ph - deg2rad(315)) < eps);
%             
%             plot(rad2deg(obj.th(iph0)),Ziplot(iph0),[LineStyle,'k'],'lineWidth',lw), grid on, hold on
%             plot(rad2deg(obj.th(iph45)),Ziplot(iph45),[LineStyle,'b'],'lineWidth',lw)
%             plot(rad2deg(obj.th(iph90)),Ziplot(iph90),[LineStyle,'r'],'lineWidth',lw)
%             plot(rad2deg(obj.th(iph135)),Ziplot(iph135),[LineStyle,'g'],'lineWidth',lw)
%             plot(-rad2deg(obj.th(iph180)),Ziplot(iph180),[LineStyle,'k'],'lineWidth',lw)
%             plot(-rad2deg(obj.th(iph225)),Ziplot(iph225),[LineStyle,'b'],'lineWidth',lw)
%             plot(-rad2deg(obj.th(iph270)),Ziplot(iph270),[LineStyle,'r'],'lineWidth',lw)
%             plot(-rad2deg(obj.th(iph315)),Ziplot(iph315),[LineStyle,'g'],'lineWidth',lw)
%             lg = legend('\phi = 0^\circ','\phi = 45^\circ','\phi = 90^\circ','\phi = 135^\circ');
%             xlabel('\theta (deg)')
%             ylabel([outputType,'(', compName, ') (',unit,')'])
%             %                 lg.String = lg.String(1:4); % Just keep the first four legend entries
%         else
%             switch cutConstant
%                 case 'x'
%                     plot(yiplot,Ziplot), grid on
%                     xlabel(yname)
%                     cutName = obj.xname;
%                 case 'y'
%                     plot(xiplot,Ziplot), grid on
%                     xlabel(xname)
%                     cutName = obj.yname;
%             end
%             ylabel([outputType,'(', compName, ') (',unit,')'])
%             
%         end
%         % Handle dynamic range here
%         if strcmp(outputType,'mag')
%             maxVal = max(Zplot(:));
%             switch scaleMag
%                 case 'dB'
%                     if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
%                         ylim([0,dynamicRange_dB] - norm.*dynamicRange_dB);
%                     else
%                         ylim([maxVal-dynamicRange_dB,maxVal]);
%                     end
%                 case 'lin'
%                     linHandle = str2func(['lin',num2str(dBscale)]);
%                     if ~norm
%                         dr = linHandle(dynamicRange_dB);
%                         if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
%                             ylim([0,dr]);
%                         else
%                             ylim([maxVal/dr,maxVal]);
%                         end
%                     end
%             end
%         end
%         if ~strcmp(obj.gridType,'DirCos')
%             cutValue = rad2deg(cutValue);
%         end
%         title([obj.coorSys, ', ',obj.polType, ' polarisation; Freq = ',num2str(freqPlot),' ', freqUnit,'; ',cutName, ' = ',num2str(cutValue), ' ',axisUnit])
%     case 'polar'
%         if ~strcmp(obj.gridType,'DirCos')
%             switch cutConstant
%                 case 'x'
%                     polarplot(deg2rad(yiplot),Ziplot), grid on
%                     cutName = obj.xname;
%                 case 'y'
%                     polarplot(deg2rad(xiplot),Ziplot), grid on
%                     cutName = obj.yname;
%             end
%             % Handle dynamic range here
%             if strcmp(outputType,'mag')
%                 maxVal = max(Zplot(:));
%                 switch scaleMag
%                     case 'dB'
%                         if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
%                             rlim([0,dynamicRange_dB] - norm.*dynamicRange_dB);
%                         else
%                             rlim([maxVal-dynamicRange_dB,maxVal]);
%                         end
%                     case 'lin'
%                         linHandle = str2func(['lin',num2str(dBscale)]);
%                         if ~norm
%                             dr = linHandle(dynamicRange_dB);
%                             if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
%                                 rlim([0,dr]);
%                             else
%                                 rlim([maxVal/dr,maxVal]);
%                             end
%                         end
%                 end
%             end
%             % Sort ou the axis
%             ax = gca;
%             ax.ThetaZeroLocation = 'top';
%             ax.ThetaDir = 'clockwise';
%         else
%             error(['Polar plots not supported for DirCos gridType'])
%         end
        
end













% %% Old code under here...
% switch obj.gridType
%     case 'UV'
%         x = obj.x;
%         y = obj.y;
%         xname = [obj.xname, ' = sin(\theta)cos(\phi)'];
%         yname = [obj.yname, ' = sin(\theta)sin(\phi)'];
%     otherwise
%         x = rad2deg(obj.x);
%         y = rad2deg(obj.y);
%         xname = [obj.xname, ' (deg)'];
%         yname = [obj.yname, ' (deg)'];
% end
% 
% 
% %% Extract plot output type
% 
% switch output
%     case {'Directivity','Gain','AxialRatio','AxialRatioInv','CO/XP','XP/CO'}
%         dr = lin10(-dynamicRange_dB);
%         if strcmp(output,'Directivity')
%             Zmat = getDirectivity(obj);
%             compName = 'Directivity';
%         elseif strcmp(output,'Gain')
%             Zmat = getGain(obj);
%             compName = 'Gain (IEEE)';
%         elseif strcmp(output,'AxialRatio')
%             [Zmat,~] = getAxialRatio(obj);
%             compName = 'Axial Ratio';
%         elseif strcmp(output,'AxialRatioInv')
%             [~,Zmat] = getAxialRatio(obj);
%             compName = 'Inverse Axial Ratio';
%         elseif strcmp(output,'CO/XP') || strcmp(output,'XP/CO')
%             switch output
%                 case 'CO/XP'
%                     num = abs(obj.E2);
%                     den = abs(obj.E1);
%                     compName = 'Co-pol/Cross-pol';
%                 case 'XP/CO'
%                     num = abs(obj.E1);
%                     den = abs(obj.E2);
%                     compName = 'Cross-pol/Co-pol';
%             end
%             % Sort out dynamic range here to avoid divide by zero
%             maxVal = max(max(num));
%             minVal = lin20(-dynamicRange_dB).*maxVal;
%             den(den < minVal) = minVal;
%             Zmat = (num./den).^2;
%         end
%         if norm, Zmat = bsxfun(@times,Zmat,1./(max(Zmat))); end
%         if strcmp(scaleMag,'dB')
%             Zmat = dB10(Zmat); unit = 'dB';
%             if norm, unit = 'dBi'; end
%         else
%             unit = 'pu';
%         end
%     case {'E1','E2'}
%         dr = lin20(-dynamicRange_dB);
%         if strcmp(output,'E1')
%             [Zmat,~,~] = getEfield(obj);
%             compName = obj.E1name;
%         elseif strcmp(output,'E2')
%             [~,Zmat,~] = getEfield(obj);
%             compName = obj.E2name;
%         end
%         if strcmp(outputType,'mag')
%             Zmat = abs(Zmat);
%             if norm, Zmat = bsxfun(@times,Zmat,1./(max(Zmat))); end
%             if strcmp(scaleMag,'dB')
%                 Zmat = dB20(Zmat); unit = 'dBV/m';
%                 if norm, unit = 'dB'; end
%             else
%                 if norm, unit = ''; else, unit = 'V/m'; end
%             end
%         elseif strcmp(outputType,'phase')
%             Zmat = angle(Zmat);
%             if norm, Zmat = bsxfun(@plus,Zmat,1./(min(Zmat))); end
%             unit = 'rad';
%             if strcmp(scalePhase,'deg'), Zmat = rad2deg(Zmat); unit = 'deg'; end
%         end
% end
% 
% % Z = Zmat(:,freqIndex);
% % freqVect = obj.freq(freqIndex);
% 
% %% Make the plots
% switch freqUnit
%     case 'Hz'
%         freqMult = 1;
%     case 'kHz'
%         freqMult = 1e-3;
%     case 'MHz'
%         freqMult = 1e-6;
%     case 'GHz'
%         freqMult = 1e-9;
%     case 'THz'
%         freqMult = 1e-12;
% end
% freqPlot = obj.freqHz(freqIndex)*freqMult;
% Zplot = Zmat(:,freqIndex);
% % Fix Dynamic Range if required
% if strcmp(outputType,'mag')
%     if strcmp(scaleMag,'dB')
%         Zplot(Zplot < (max(Zplot) - dynamicRange_dB)) = max(Zplot) - dynamicRange_dB;
%     else
%         Zplot(Zplot < max(Zplot)*dr) = max(Zplot)*dr;
%     end
% end
% 
% % keyboard;
% 
% if isempty(step)
%     xMat = reshape(x,obj.Ny,obj.Nx);
%     yMat = reshape(y,obj.Ny,obj.Nx);
%     X = xMat;
%     Y = yMat;
% elseif ~strcmp(plotType,'3D') % && ~isempty(cutValue)
%     [xySort,iSort] = unique([x,y],'rows');
%     x = xySort(:,1);
%     y = xySort(:,2);
%     Zplot = Zplot(iSort);
%     xvect = min(x):step:max(x);
%     yvect = min(y):step:max(y);
%     [X,Y] = meshgrid(xvect,yvect);
%     Zplot = griddata(x,y,Zplot,X,Y,'linear');
%     x = X(:);
%     y = Y(:);
% end
% 
% %     figure
% switch plotType
%     case '3D'
%         % Use the MATLAB antennas toolbox plotting function
%         patternCustom(Zplot,rad2deg(obj.th),rad2deg(obj.ph));
%         title([obj.coorSys, ', ',obj.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,'); Freq = ',num2str(freqPlot),' ', freqUnit])
%     case '2D'
%         surf(X,Y,reshape(Zplot,size(X)),'EdgeColor','Interp','FaceColor','Interp')
%         xlabel(xname)
%         ylabel(yname)
%         view([0,90])
%         axis equal
%         xlim([min(x),max(x)])
%         ylim([min(y), max(y)])
%         colorbar
%         title([obj.coorSys, ', ',obj.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,'); Freq = ',num2str(freqPlot),' ', freqUnit])
%     case 'cartesian'
%         % ToDo
%         lw = LineWidth;   
%         if isempty(cutValue)
%             % Find the principle cuts
%             iph0 = find(abs(obj.ph - 0) < eps);
%             iph45 = find(abs(obj.ph - deg2rad(45)) < eps);
%             iph90 = find(abs(obj.ph - deg2rad(90)) < eps);
%             iph135 = find(abs(obj.ph - deg2rad(135)) < eps);
%             iph180 = find(abs(obj.ph - deg2rad(180)) < eps);
%             iph225 = find(abs(obj.ph - deg2rad(225)) < eps);
%             iph270 = find(abs(obj.ph - deg2rad(270)) < eps);
%             iph315 = find(abs(obj.ph - deg2rad(315)) < eps);
%             
%             plot(rad2deg(obj.th(iph0)),Zplot(iph0),[LineStyle,'k'],'lineWidth',lw), grid on, hold on
%             plot(rad2deg(obj.th(iph45)),Zplot(iph45),[LineStyle,'b'],'lineWidth',lw)
%             plot(rad2deg(obj.th(iph90)),Zplot(iph90),[LineStyle,'r'],'lineWidth',lw)
%             plot(rad2deg(obj.th(iph135)),Zplot(iph135),[LineStyle,'g'],'lineWidth',lw)
%             plot(-rad2deg(obj.th(iph180)),Zplot(iph180),[LineStyle,'k'],'lineWidth',lw)
%             plot(-rad2deg(obj.th(iph225)),Zplot(iph225),[LineStyle,'b'],'lineWidth',lw)
%             plot(-rad2deg(obj.th(iph270)),Zplot(iph270),[LineStyle,'r'],'lineWidth',lw)
%             plot(-rad2deg(obj.th(iph315)),Zplot(iph315),[LineStyle,'g'],'lineWidth',lw)
%             lg = legend('\phi = 0^\circ','\phi = 45^\circ','\phi = 90^\circ','\phi = 135^\circ');
%             xlabel('\theta (deg)')
%             ylabel([outputType,'(', compName, ') (',unit,')'])
%             %                 lg.String = lg.String(1:4); % Just keep the first four legend entries
%         else
%             if cutConstant == 'x'
%                 iCut = find(abs(x - cutValue) < eps);
%                 plot(y(iCut),Zplot(iCut),LineStyle,'lineWidth',lw), grid on, hold on
%                 xlabel(yname)
%             elseif cutConstant == 'y'
%                 iCut = find(abs(y - cutValue) < eps);
%                 plot(x(iCut),Zplot(iCut),LineStyle,'lineWidth',lw), grid on, hold on
%                 xlabel(xname)
%             end
%             ylabel([outputType,'(', compName, ') (',unit,')'])
%             
%         end
%         title([obj.coorSys, ', ',obj.polType, ' polarisation; Freq = ',num2str(freqPlot),' ', freqUnit])
%     case 'polar'
%         % ToDo
% 
% end


% end

end