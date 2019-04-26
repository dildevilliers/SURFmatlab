function [Zi] = interpolateGrid(obj,output,xi,yi,varargin)

%% Check inputs
narginchk(4,6);

parseobj = inputParser;
parseobj.FunctionName = 'interpolateGrid';

typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'interpolateGrid','obj',1);
addRequired(parseobj,'obj',typeValidationObj);

expectedoutput = {'Directivity','Gain','E1','E2','E3','AxialRatio','AxialRatioInv','CO_XP','XP_CO','U','W'};
addRequired(parseobj,'output', @(x) any(validatestring(x,expectedoutput)));

typeValidationXY = @(x) validateattributes(x,{'numeric'},{'real','nonempty'},'interpolateGrid');
addRequired(parseobj,'xi',typeValidationXY);
addRequired(parseobj,'yi',typeValidationXY);

% expectedgridType = {'PhTh','DirCos','AzEl','ElAz','TrueView','ArcSin'};
% addOptional(parseobj,'gridType',obj.gridType, @(x) any(validatestring(x,expectedgridType)));

typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'interpolateGrid','freqIndex');
addOptional(parseobj,'freqIndex',1,typeValidationFreq);
% addParameter(parseobj,'freqIndex',1,typeValidationFreq);

expectedhemisphere = {'top','bot'};
addOptional(parseobj,'hemisphere','top', @(x) any(validatestring(x,expectedhemisphere)));
% addParameter(parseobj,'hemisphere','top', @(x) any(validatestring(x,expectedhemisphere)));

parse(parseobj,obj,output,xi,yi,varargin{:});

% gridType = parseobj.Results.gridType;
freqIndex = parseobj.Results.freqIndex;
hemisphere = parseobj.Results.hemisphere;


%% Main code

gridType = obj.gridType;

% Evaluate the field on the base grid - this is where the output function
% should be best suited for interpolation
obj = obj.grid2Base;
% Shift to -180:180 range (if applicable) - this is where the DirCos spits
% everything out after transforming
if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz') %|| strcmp(obj.gridType,'AzAlt')
    obj = obj.setXrange('sym');
end
% Get xi and yi in the base gridType, and on the [-180,180] x-domain for the
% angular grids
grid2DirCoshandle = str2func([gridType,'2DirCos']);
[ui,vi,wi] = grid2DirCoshandle(xi,yi);
% Check for bottom hemisphere plot - fix the w to be the negative root
if (strcmp(gridType,'DirCos') || strcmp(gridType,'ArcSin')) && strcmp(hemisphere,'bot')
    wi = -wi;
end
DirCos2baseHandle = str2func(['DirCos2',obj.gridType]);
[xi_bGT,yi_bGT] = DirCos2baseHandle(ui,vi,wi);
% Find the invalid points included by the external meshgrid 
valAngi = sqrt(ui.^2 + vi.^2) <= 1;
% Sort out the TrueView special case invalid points
if strcmp(gridType,'TrueView')
    valAngi = sqrt((xi./pi).^2 + (yi./pi).^2) <= 1;
end

% Get the valid angle positions - already in baseGrid here, but shifted to
% +- 180 degrees
% Get the indexes only, no later reshaping done, different from the grids
% required for plotting in plot.m
if strcmp(gridType,'DirCos') || strcmp(gridType,'ArcSin')  
    grid2DirCoshandleBase = str2func([obj.gridType,'2DirCos']);
    [~,~,w] = grid2DirCoshandleBase(obj.x,obj.y);
    if strcmp(hemisphere,'top')
        valAng = find(w >= 0);
    elseif strcmp(hemisphere,'bot')
        valAng = find(w <= 0);
    end
else
    valAng = 1:obj.Nang;
end
valAng = valAng(:);

% Extract the outputs on the base grid
if strcmp(output,'E1')
    [Zfreq,~,~] = getEfield(obj);
elseif strcmp(output,'E2')
    [~,Zfreq,~] = getEfield(obj);
elseif strcmp(output,'E3')
    [~,~,Zfreq] = getEfield(obj);
else
    outputHandle = str2func(['get',output]);
    Zfreq = outputHandle(obj);
end

% Select frequency of interest
Z = Zfreq(:,freqIndex);

xVal = obj.x(valAng);
yVal = obj.y(valAng);
zVal = Z(valAng);

edgeAngExtent_deg = 16;
% Extend grid past -180 and +180 for interpolation across the axis
if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
    tol = deg2rad(2); % Check for points close to the edges in azimuth
    if abs(min(xVal) + pi) < tol && abs(max(xVal) - pi) < tol 
        edgeAngDeg = 180 - edgeAngExtent_deg;
        iNeg = find(xVal > deg2rad(edgeAngDeg));
        iPos = find(xVal < deg2rad(-edgeAngDeg));
        xVal = [xVal(iNeg)-2*pi;xVal;xVal(iPos)+2*pi];
        yVal = [yVal(iNeg);yVal;yVal(iPos)];
        zVal = [zVal(iNeg);zVal;zVal(iPos)];
    end
%     % Also extend the y-axis
%     if abs(min(yVal) - 0) < tol
%         edgeAngDeg = edgeAngExtent_deg;
%         iNeg = find(xVal < 0 & yVal < deg2rad(edgeAngDeg));
%         iPos = find(xVal > 0 & yVal < deg2rad(edgeAngDeg));
%         xVal = [xVal;xVal(iNeg)+pi;xVal(iPos)-pi];
%         yVal = [yVal;-yVal(iNeg);-yVal(iPos)];
%         zVal = [zVal;zVal(iNeg);zVal(iPos)];
%     end
%     if abs(max(yVal) - pi) < tol
%         edgeAngDeg = 180 - edgeAngExtent_deg;
%         iNeg = find(xVal < 0 & yVal > deg2rad(edgeAngDeg));
%         iPos = find(xVal > 0 & yVal > deg2rad(edgeAngDeg));
%         xVal = [xVal;xVal(iNeg)+pi;xVal(iPos)-pi];
%         yVal = [yVal;2*pi-yVal(iNeg);2*pi-yVal(iPos)];
%         zVal = [zVal;zVal(iNeg);zVal(iPos)];
%     end
    
end

% Remove duplicate differing values completely from the set - interpolate
% over them.  This happens at poles for certain coordinate projections
% (like at th = 180 in Ludwig3 for instance, of th=0 and 180 for spherical)
% First find duplicate domain values
[~,iUnique] = unique([xVal,yVal],'rows');
removePoints = [];
if length(iUnique) < length(xVal)
    iRepeated = setdiff((1:length(xVal)).',iUnique);
    repeatedSet = [xVal(iRepeated),yVal(iRepeated),zVal(iRepeated)];
    % Find the repeated values
    repVals = unique(repeatedSet(:,1:2),'rows');
    % Test if different z-values occur
    for ii = 1:length(repVals(:,1))
        currentValRowIndex = find(ismember(repeatedSet(:,1:2),repVals(ii,:),'rows'));
        currentZ = repeatedSet(currentValRowIndex,3);
        if ~all(currentZ == currentZ(1))
            removePoints = [removePoints;repeatedSet(currentValRowIndex,1:2)];
        end
    end
end
if numel(removePoints) > 0
    iRemove = find(ismember([xVal,yVal],removePoints,'rows'));
    xVal(iRemove) = [];
    yVal(iRemove) = [];
    zVal(iRemove) = [];
end

% Build the interpolant on the base grid at the valid angles
if obj.isGridUniform
try
    NyVal = length(unique(yVal));
    NxVal = length(unique(xVal));
    XVal = reshape(xVal,NyVal,NxVal);
    YVal = reshape(yVal,NyVal,NxVal);
    ZVal = reshape(zVal,NyVal,NxVal);
    Zf = griddedInterpolant(XVal',YVal',ZVal.','linear');
catch ME
    % Grid did not work... Go to scatter
end
end
if ~exist('Zf','var')
    warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
    Zf = scatteredInterpolant(xVal,yVal,zVal,'linear');
    warning('on','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
end
% Get the values on the scattered set of inputs
Zi = Zf(xi_bGT,yi_bGT);
Zi(~valAngi) = NaN;
end
