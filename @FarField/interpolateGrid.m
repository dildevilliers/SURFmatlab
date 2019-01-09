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
if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
    obj = obj.setXrange('sym');
end
% Get xi and yi in the base gridType, and on the [-180,180] x-domain for the
% angular grids
grid2DirCoshandle = str2func(['FarField.',gridType,'2DirCos']);
[ui,vi,wi] = grid2DirCoshandle(xi,yi);
% Check for bottom hemisphere plot - fix the w to be the negative root
if (strcmp(gridType,'DirCos') || strcmp(gridType,'ArcSin')) && strcmp(hemisphere,'bot')
    wi = -wi;
end
DirCos2baseHandle = str2func(['FarField.DirCos2',obj.gridType]);
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
    grid2DirCoshandleBase = str2func(['FarField.',obj.gridType,'2DirCos']);
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

% Build the interpolant on the base grid at the valid angles
Zf = scatteredInterpolant(obj.x(valAng),obj.y(valAng),Z(valAng),'linear');
% Get the values on the scattered set of inputs
Zi = Zf(xi_bGT,yi_bGT);
Zi(~valAngi) = NaN;
end
