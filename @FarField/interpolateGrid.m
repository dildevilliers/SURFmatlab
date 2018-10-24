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

expectedgridType = {'PhTh','UV','AzEl','ElAz','TrueView','ArcSin'};
addOptional(parseobj,'gridType',obj.gridType, @(x) any(validatestring(x,expectedgridType)));

typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'interpolateGrid','freqIndex');
addOptional(parseobj,'freqIndex',1,typeValidationFreq);

parse(parseobj,obj,output,xi,yi,varargin{:});

gridType = parseobj.Results.gridType;
freqIndex = parseobj.Results.freqIndex;

%% Main code

obj = obj.resetToBase;

% Get xi and yi in the base gridType
if ~strcmp(obj.gridType,gridType)
    grid2UVhandle = str2func(['FarField.',gridType,'2UV']);
    [u,v,w] = grid2UVhandle(xi,yi);
    
    UV2baseHandle = str2func(['FarField.UV2',obj.gridTypeBase]);
    [xi_bGT,yi_bGT] = UV2baseHandle(u,v,w);
else
    xi_bGT = xi;
    yi_bGT = yi;
end

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

% Build the interpolant on the base grid
Zf = scatteredInterpolant(obj.x,obj.y,Z,'linear','linear');
% Get the values on the scattered set of inputs
Zi = Zf(xi_bGT,yi_bGT);

end
