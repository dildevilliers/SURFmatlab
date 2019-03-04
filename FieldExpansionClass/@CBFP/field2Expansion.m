function CBFPobj = field2Expansion(FFobj,tol,Nbasis,i_FM)

if nargin < 2
    tol = 1e-100;
    Nbasis = inf;  
    i_FM = [];
elseif nargin < 3
    Nbasis = inf;
    i_FM = [];
elseif nargin < 4
    i_FM = [];
end

%% Extract useful properties

Nf = FFobj.Nf;
E1 = FFobj.E1;
E2 = FFobj.E2;
Nang = FFobj.Nang;

%% Basic error checking

% if FF(1).Nf > 1 assert(all(ismember(reshape([FF.freq],length(FF(1).freq),length(FF))',FF(1).freq,'rows')),'All FF.freq fields should be identical'); end
% assert(all([FF.Nth] == FF(1).Nth),'All FF.Nth fields should be identical')
% assert(all([FF.Nph] == FF(1).Nph),'All FF.Nph fields should be identical')

%% Lengths

% Nx = length(FF); % parametric variation
% Nf = length(FF(1).freq);
% Np = length(FF(1).th(:,1));

%% Determine FF indices from which CBFP matrix is built

if isfield(i_FM,'f') 
    range_f = reshape(i_FM.f,1,length(i_FM.f)); 
else
    range_f = 1:Nf; 
end
    
%% Build CBFP matrices

% Get the CBFP matrix
FM_E1 = [];
FM_E2 = [];
FM_E1 = [FM_E1 E1(:,range_f)];
FM_E2 = [FM_E2 E2(:,range_f)];

FM = [FM_E1;FM_E2];

% Get the SVD of the FM matrix
[U,S,V] = svd(FM,'econ');

Snorm = S(1,1);
sigma_n = diag(S)./Snorm;   % Normalization factor for sigma - this leaves an option to zero one of the field components if required and get rid of noise

% Sort out the SVD matrix sizes due to reduced number of significant singular values
if Nbasis < inf
    NR = Nbasis;
else
    NR = length(find(sigma_n > tol));
end

if NR > 0
    UR = U(:,1:NR);
    VR = V(:,1:NR);
    R = bsxfun(@rdivide,FM*VR,(sigma_n(1:NR)'*Snorm));
end

% Calculate CBFP weights for all requested input points (not just at the CBFP support points)
W = pinv(R)*FM;

%% Build the Object

% Handle basis functions
x = FFobj.x;
y = FFobj.y;
freq = 1;
Prad = ones(size(freq)).*4*pi;
radEff = ones(size(freq));
coorSys = FFobj.coorSys;
polType = FFobj.polType;
gridType = FFobj.gridType;
freqUnit = FFobj.freqUnit;

for ii = 1:NR
    
    basis_E1 = R(1:Nang,ii);
    basis_E2 = R(Nang+1:end,ii);
    basis_E3 = ones([Nang,1]);
    
    FFbasis = FarField(x,y,basis_E1,basis_E2,basis_E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
    
    basis(1,ii) = FFbasis;
end

% FF_basis = FarField.setEnames(FF_basis);
% FF_basis = FarField.setXYnames(FF_basis);
% FF_basis = FarField.setPhTh(FF_basis);
% FF_basis = FarField.setFreq(FF_basis);
% FF_basis = FarField.setBase(FF_basis);

% CBFP Object

coeffs = W;
nBasis = NR;
nCoeffs = length(W(:,1));

CBFPobj = CBFP(nBasis,nCoeffs,UR,S,VR,basis,coeffs,sigma_n);
end
