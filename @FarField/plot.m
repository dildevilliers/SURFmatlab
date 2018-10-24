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
%   'CO_XP' | 'XP_CO'
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

typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'plot','obj',1);
addRequired(parseobj,'obj',typeValidationObj);

typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'plot','freqIndex');
addParameter(parseobj,'freqIndex',1,typeValidationFreq);

typeValidationnorm = @(x) validateattributes(x,{'numeric'},{'binary','nonempty','numel',1},'plot','norm');
addParameter(parseobj,'norm',false,typeValidationnorm );

typeValidationDR = @(x) validateattributes(x,{'numeric'},{'real','positive','nonempty','numel',1},'plot','dynamicRange_dB');
addParameter(parseobj,'dynamicRange_dB',40,typeValidationDR );

expectedplotType = {'3D','2D','polar','cartesian'};
addParameter(parseobj,'plotType','3D', @(x) any(validatestring(x,expectedplotType)));

expectedoutput = {'Directivity','Gain','E1','E2','E3','AxialRatio','AxialRatioInv','CO_XP','XP_CO','W','U'};
addParameter(parseobj,'output','Directivity', @(x) any(validatestring(x,expectedoutput)));

expectedoutputType = {'mag','phase'};
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

addParameter(parseobj,'LineStyle','-');

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

%% Sort out the plot grid and names

if isempty(step)
    xi = obj.x;
    yi = obj.y;
    NxPlot = obj.Nx;
    NyPlot = obj.Ny;
    X = reshape(xi,NyPlot,NxPlot);
    Y = reshape(yi,NyPlot,NxPlot);
    
    % Extract the outputs on the base grid
    if strcmp(output,'E1')
        [Zi,~,~] = getEfield(obj);
    elseif strcmp(output,'E2')
        [~,Zi,~] = getEfield(obj);
    elseif strcmp(output,'E3')
        [~,~,Zi] = getEfield(obj);
    else
        outputHandle = str2func(['get',output]);
        Zi = outputHandle(obj);
    end
    Zi = Zi(:,freqIndex);

else
    % Get the plot points from the step information
    if strcmp(obj.gridType,'UV') || strcmp(obj.gridType,'ArcSin')
        step = sind(step);
    else
        step = deg2rad(step);
    end
    xvect = min(obj.x):step:max(obj.x);
    yvect = min(obj.y):step:max(obj.y);
    [X,Y] = meshgrid(xvect,yvect);
    NxPlot = numel(xvect);
    NyPlot = numel(yvect);
    switch plotType
        case {'3D','2D'}
            xi = X(:);
            yi = Y(:);
        case {'cartesian','polar'}
            switch cutConstant
                case 'x'
                    xi = Y(:);
                case 'y'
                    xi = X(:);
            end
    end
    [Zi] = interpolateGrid(obj,output,xi,yi,obj.gridType,freqIndex);
end

% Assign axis names
switch obj.gridType
    case 'UV'
        xplot = xi;
        yplot = yi;
        xname = [obj.xname, ' = sin(\theta)cos(\phi)'];
        yname = [obj.yname, ' = sin(\theta)sin(\phi)'];
    otherwise
%         xvect = rad2deg(xvect);
%         yvect = rad2deg(yvect);
        X = rad2deg(X);
        Y = rad2deg(Y);
        xplot = rad2deg(xi);
        yplot = rad2deg(yi);
        xname = [obj.xname, ' (deg)'];
        yname = [obj.yname, ' (deg)'];
end

%% Condition outputs
% Phase results
if (strcmp(output,'E1') || strcmp(output,'E2') || strcmp(output,'E3')) && strcmp(outputType,'phase')
    Zplot = angle(Zi);
    if norm, Zplot = Zplot - max(Zplot); end
    unit = 'V/m (rad)';
    if strcmp(scalePhase,'deg')
        Zplot = rad2deg(Zplot); 
        unit = 'V/m (deg)';
    end
else    % Power results
    Zplot = abs(Zi);
    switch output
        case {'Directivity','Gain','AxialRatio','AxialRatioInv','CO_XP','XP_CO','W','U'}
            dBscale = 10;
            unit = '';
            compName = output;
        case {'E1','E2','E3'}
            dBscale = 20;
            unit = 'V/m ';
            compName = obj.([output,'name']);
    end
    linHandle = str2func(['lin',num2str(dBscale)]);
    dr = linHandle(dynamicRange_dB);
    if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
        minVal = min(Zplot);
        maxVal = dr*minVal;
        Zplot(Zplot > maxVal) = maxVal;
    else
        maxVal = max(Zplot);
        minVal = maxVal/dr;
        Zplot(Zplot < minVal) = minVal;
    end
    if norm, Zplot = Zplot./maxVal; end
    if strcmp(scaleMag,'dB')
        dBHandle = str2func(['dB',num2str(dBscale)]);
        Zplot = dBHandle(Zplot); 
        unit = [unit, '(dB)'];
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
        % Use the MATLAB antennas toolbox plotting function
        
        % Need the th and ph angles for this so get them
        if strcmp(obj.gridType,obj.gridTypeBase)
            phplot = xplot;
            thplot = yplot;
        else
            grid2UVhandle = str2func(['FarField.',obj.gridType,'2UV']);
            [u,v,w] = grid2UVhandle(xi,yi);
            [phi,thi] = FarField.UV2PhTh(u,v,w);
            [angMat,IA] = unique([phi,thi],'rows');
            phplot = rad2deg(angMat(:,1));
            thplot = rad2deg(angMat(:,2));
            Zplot = Zplot(IA);
        end
        patternCustom(Zplot,thplot,phplot);
        title([obj.coorSys, ', ',obj.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,'); Freq = ',num2str(freqPlot),' ', freqUnit])
        
    case '2D'
%         keyboard;
%         Xplot = reshape(xplot,NxPlot,NyPlot);
%         Yplot = reshape(yplot,NxPlot,NyPlot);
%         [Zplot] = griddata(reshape(xplot,NxPlot,NyPlot),reshape(yplot,NxPlot,NyPlot),reshape(Zplot,NxPlot,NyPlot),X,Y);
%         Xplot = X;
%         Yplot = Y;
%         iInvisible = find((Xplot.^2 + Yplot.^2) > 1);
%         Zplot(iInvisible) = NaN;
%         surf(Xplot,Yplot,reshape(Zplot,size(Xplot)),'EdgeColor','Interp','FaceColor','Interp')
Xplot = X;
Yplot = Y;
Zplot = reshape(Zplot,NyPlot,NxPlot);
% keyboard;
        surf(Xplot,Yplot,Zplot,'EdgeColor','Interp','FaceColor','Interp')
%         figure
%         plot3(xplot,yplot,Zplot,'.')
        xlabel(xname)
        ylabel(yname)
        view([0,90])
        axis equal
        xlim([min(xplot),max(xplot)])
        ylim([min(yplot), max(yplot)])
        colorbar
        title([obj.coorSys, ', ',obj.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,'); Freq = ',num2str(freqPlot),' ', freqUnit])
    case 'cartesian'
        % ToDo
        lw = LineWidth;   
        if isempty(cutValue)
            % Find the principle cuts
            iph0 = find(abs(obj.ph - 0) < eps);
            iph45 = find(abs(obj.ph - deg2rad(45)) < eps);
            iph90 = find(abs(obj.ph - deg2rad(90)) < eps);
            iph135 = find(abs(obj.ph - deg2rad(135)) < eps);
            iph180 = find(abs(obj.ph - deg2rad(180)) < eps);
            iph225 = find(abs(obj.ph - deg2rad(225)) < eps);
            iph270 = find(abs(obj.ph - deg2rad(270)) < eps);
            iph315 = find(abs(obj.ph - deg2rad(315)) < eps);
            
            plot(rad2deg(obj.th(iph0)),Zplot(iph0),[LineStyle,'k'],'lineWidth',lw), grid on, hold on
            plot(rad2deg(obj.th(iph45)),Zplot(iph45),[LineStyle,'b'],'lineWidth',lw)
            plot(rad2deg(obj.th(iph90)),Zplot(iph90),[LineStyle,'r'],'lineWidth',lw)
            plot(rad2deg(obj.th(iph135)),Zplot(iph135),[LineStyle,'g'],'lineWidth',lw)
            plot(-rad2deg(obj.th(iph180)),Zplot(iph180),[LineStyle,'k'],'lineWidth',lw)
            plot(-rad2deg(obj.th(iph225)),Zplot(iph225),[LineStyle,'b'],'lineWidth',lw)
            plot(-rad2deg(obj.th(iph270)),Zplot(iph270),[LineStyle,'r'],'lineWidth',lw)
            plot(-rad2deg(obj.th(iph315)),Zplot(iph315),[LineStyle,'g'],'lineWidth',lw)
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
            
        end
        title([obj.coorSys, ', ',obj.polType, ' polarisation; Freq = ',num2str(freqPlot),' ', freqUnit])
    case 'polar'
        % ToDo

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