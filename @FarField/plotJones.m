function [] = plotJones(FF1,FF2,varargin)

% function [] = plotJones(FF1,FF2,varargin)
% plots the Jones matrix type representation of the farfields specified in
% FF1 and FF2, where it is assumed they are calculated in the same basis
% and represent the 2 polarizations
% name, value are name value pairs, and can be the following:
%
% freqIndex is the index of the frequency to be plotted (default 1)
%
% dynamicRange_dB is a (positive) dB value for the magnitude plot dynamic
% range (40)

%% Parse input
parseobj = inputParser;
parseobj.FunctionName = 'plotJones';

typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'plotJones','obj',1);
addRequired(parseobj,'FF1',typeValidationObj);
addRequired(parseobj,'FF2',typeValidationObj);

typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'plotJones','freqIndex');
addParameter(parseobj,'freqIndex',1,typeValidationFreq);

typeValidationDR = @(x) validateattributes(x,{'numeric'},{'real','positive','nonempty','numel',1},'plotJones','dynamicRange_dB');
addParameter(parseobj,'dynamicRange_dB',40,typeValidationDR );

parse(parseobj, FF1, FF2, varargin{:});

freqIndex = parseobj.Results.freqIndex;
dynamicRange_dB = parseobj.Results.dynamicRange_dB;


%% Plot the result

if ~isGridEqual(FF1,FF2)
    error('Base grids should be identical for the two input fields');
else
    subplot(2,2,1)
    plot(FF1,'output','E2','outputType','mag','plotType','2D','scaleMag','dB','norm',1,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex)
    title('J_{11}')
    subplot(2,2,2)
    plot(FF1,'output','E1','outputType','mag','plotType','2D','scaleMag','dB','norm',1,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex)
    title('J_{12}')
    subplot(2,2,3)
    plot(FF2,'output','E1','outputType','mag','plotType','2D','scaleMag','dB','norm',1,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex)
    title('J_{21}')
    subplot(2,2,4)
    plot(FF2,'output','E2','outputType','mag','plotType','2D','scaleMag','dB','norm',1,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex)
    title('J_{22}')
end