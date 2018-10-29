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

parse(parseobj,obj,output,xi,yi,varargin{:});

% gridType = parseobj.Results.gridType;
freqIndex = parseobj.Results.freqIndex;

%% Main code

gridType = obj.gridType;

% Evaluate the field on the base grid - this is where the output function
% should be best suited for interpolation
obj = obj.grid2Base;

% Get xi and yi in the base gridType
if ~strcmp(obj.gridType,gridType)
    % Shift to -180:180 range (if applicable)
    obj = obj.xRange180180;

    grid2DirCoshandle = str2func(['FarField.',gridType,'2DirCos']);
    [u,v,w] = grid2DirCoshandle(xi,yi);
    
    DirCos2baseHandle = str2func(['FarField.DirCos2',obj.gridTypeBase]);
    [xi_bGT,yi_bGT] = DirCos2baseHandle(u,v,w);
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
Zf = scatteredInterpolant(obj.x,obj.y,Z,'natural');
% Get the values on the scattered set of inputs
Zi = Zf(xi_bGT,yi_bGT);

end
