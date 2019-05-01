function [] = plotPrincipleCuts(FF,varargin)

% [] = plotPrincipleCuts(FF,varargin)
% plots the principle cuts of the farfield specified in FF
% name, value are name value pairs, and can be the following:
%
% freqIndex is the index of the frequency to be plotted (default 1)
%
% dynamicRange_dB is a (positive) dB value for the magnitude plot dynamic
% range (40)
%
% norm is a boolean (false) to normalize to maximum magnitude
%
% plotType can be:
%    ('cartesian') | 'polar'
%
% MainPolNr can be
% (2) | 1 - Indicates which component to plot in solid lines
%
% output can be:
%   ('Directivity') | 'Gain' | 'E1' | 'E2' | 'AxialRatio' | 'AxialRatioInv'
%   'CO_XP' | 'XP_CO' | 'W' | 'U'
%   (For Power quantities only one plot is generated, for the rest both will be 
%   plotted, with the first (solid line) determined by the input value)

%% Parse input
parseobj = inputParser;
parseobj.FunctionName = 'plotPrincipleCuts';

typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'plotPrincipleCuts','obj',1);
addRequired(parseobj,'FF',typeValidationObj);

typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'plotPrincipleCuts','freqIndex');
addParameter(parseobj,'freqIndex',1,typeValidationFreq);

typeValidationDR = @(x) validateattributes(x,{'numeric'},{'real','positive','nonempty','numel',1},'plotPrincipleCuts','dynamicRange_dB');
addParameter(parseobj,'dynamicRange_dB',40,typeValidationDR );

typeValidationnorm = @(x) validateattributes(x,{'logical','numeric'},{'binary','nonempty','numel',1},'plot','norm');
addParameter(parseobj,'norm',false,typeValidationnorm );

expectedplotType = {'polar','cartesian'};
addParameter(parseobj,'plotType','cartesian', @(x) any(validatestring(x,expectedplotType)));

expectedoutput = {'Directivity','Gain','E1','E2','E3','AxialRatio','AxialRatioInv','CO_XP','XP_CO','W','U'};
addParameter(parseobj,'output','Directivity', @(x) any(validatestring(x,expectedoutput)));

parse(parseobj, FF, varargin{:});

freqIndex = parseobj.Results.freqIndex;
dynamicRange_dB = parseobj.Results.dynamicRange_dB;
norm = parseobj.Results.norm;
plotType = parseobj.Results.plotType;
output = parseobj.Results.output;


%% Plot the result

% Estimate a nice step size
step = median(diff(unique(FF.thBase)));
if strcmp(FF.gridType,'DirCos') || strcmp(FF.gridType,'ArcSin')
    step = asin(step);
end
step = rad2deg(step);

% Output control
plotSec = true;
switch output
    case {'Directivity','Gain','W','U'}
        plotSec = false;
        Emain = output;
        ylabText = ['|',output,'| (dB)'];
    case 'E1'
        Emain = 'E1';
        Esec = 'E2';
        ylabText = ['|',FF.E1name,'| (-); |',FF.E2name,'| (--) (dB)' ];
    case 'E2'
        Emain = 'E2';
        Esec = 'E1';
        ylabText = ['|',FF.E2name,'| (-); |',FF.E1name,'| (--) (dB)' ];
    case 'AxialRatio'
        Emain = 'AxialRatio';
        plotSec = false;
        ylabText = ['|AR| (dB)' ];
    case 'AxialRatioInv'
        Emain = 'AxialRatioInv';
        plotSec = false;
        ylabText = ['|AR| (dB)' ];
    case 'CO_XP'
        Emain = 'CO_XP';
        plotSec = false;
        ylabText = ['|CO/XP| (dB)' ];
    case 'XP_CO'
        Emain = 'XP_CO';
        plotSec = false;
        ylabText = ['|XP/CO| (dB)' ];
end

figure
switch FF.gridType
    case{'PhTh','AzEl','ElAz'}
        % Shift the pattern onto a symmetrical grid
        if ~FF.isGridUniform
            FF = currentForm2Base(FF,step);
        else
            FF = FF.reset2Base;
        end
        FF = FF.setXrange('sym');
        FF = FF.setYrange(360);
        xVal1 = 0;
        xVal2 = 90;
        xVal3 = 45;
        
        % Main component
        plot(FF,'output',Emain,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
            'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal1),...
            'LineStyle','-','Color','k')
        plot(FF,'output',Emain,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
            'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal2),...
            'LineStyle','-','Color','r')
        plot(FF,'output',Emain,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
            'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal3),...
            'LineStyle','-','Color','b')
        
        if plotSec
            % 2nd Component
            plot(FF,'output',Esec,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal1),...
                'LineStyle','--','Color','k')
            hold on
            plot(FF,'output',Esec,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal2),...
                'LineStyle','--','Color','r')
            plot(FF,'output',Esec,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal3),...
                'LineStyle','--','Color','b')
        end
        % Add the legend
        xUnit = '^\circ';
        legend([FF.xname,'=',num2str(xVal1),xUnit],[FF.xname,'=',num2str(xVal2),xUnit],[FF.xname,'=',num2str(xVal3),xUnit])

    otherwise
        xVal = 0;
        yVal = 0;
        % Main component
        plot(FF,'output',Emain,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
            'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal),...
            'cutConstant','x','LineStyle','-','Color','k')
        plot(FF,'output',Emain,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
            'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(yVal),...
            'cutConstant','y','LineStyle','-','Color','r')
        if plotSec
            % 2nd Component
            plot(FF,'output',Esec,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal),...
                'cutConstant','x','LineStyle','--','Color','k')
            plot(FF,'output',Esec,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(yVal),...
                'cutConstant','y','LineStyle','--','Color','r')
        end

        % Add the legend
        legend([FF.xname,'=',num2str(xVal)],[FF.yname,'=',num2str(yVal)])
end
if isequal(plotType,'cartesian')
    ylabel(ylabText)
end

% Remove the cut value from the title
h = gca;
titTextFull = h.Title.String;
titTextCell = strsplit(titTextFull,'Hz');
title([titTextCell{1},'Hz'])






















