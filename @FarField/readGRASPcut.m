function [FF] = readGRASPcut(pathName,nr_freq,nr_cuts)

% [FF] = readGRASPcut(pathName)
% Loads a GRASP generated farfield source file pathName.cut.
%
%
% Inputs:
% path_name - Full path and filename string
% nr_freq - number of frequency points
% nr_cuts - number of cuts taken
% Outputs:
% FF - standard farfield object
%
% Dirk de Villiers
% Created: 2019-04-22


%Open the data file
global fid;
[fid, message] = fopen([pathName,'.cut'], 'rt');
if (fid==-1)
   error(['Unable to open data file ' fileName '!']);
end


%===================================================================
% Load data for pre-allocation
%===================================================================
% Skip over text line
form= '%*s %*s %*s %*s';
dummy = textscan(fid, form, 1);

% Read info line
form= '%f %f %d %f %d %d %d';
cut_info = textscan(fid, form, 1);
V_INI = cut_info{1};
V_INC = cut_info{2};
V_NUM = cut_info{3};
C = cut_info{4};
ICOMP = cut_info{5};
ICUT = cut_info{6};
% NCOMP = cut_info{7};

% Preallocate
[th_deg,ph_deg] = deal(zeros(V_NUM*nr_cuts,1));
[E1,E2] = deal(zeros(V_NUM*nr_cuts,nr_freq));

for ff = 1:nr_freq
    for cc = 1:nr_cuts
        if ff == 1 % Only do once
            x_1cut = ones(V_NUM,1).*C;
            y_1cut = (V_INI:V_INC:(V_INC*(V_NUM - 1) + V_INI)).';
            
            if ICUT == 1
                ph_deg(((cc-1)*V_NUM + 1):cc*V_NUM) = x_1cut;
                th_deg(((cc-1)*V_NUM + 1):cc*V_NUM) = y_1cut;
            elseif ICUT == 2
                th_deg(((cc-1)*V_NUM + 1):cc*V_NUM) = x_1cut;
                ph_deg(((cc-1)*V_NUM + 1):cc*V_NUM) = y_1cut;
            end
        end
        % Read cut data
        form= '%f %f %f %f';
        cut_data = textscan(fid, form, V_NUM);
        E1(((cc-1)*V_NUM + 1):cc*V_NUM,ff) = cut_data{1} + 1i.*cut_data{2};
        E2(((cc-1)*V_NUM + 1):cc*V_NUM,ff) = cut_data{3} + 1i.*cut_data{4};
        
        form = '%*s %*s %*s %*s';
        dummy = textscan(fid, form, 1);
        form= '%f %f %d %f %d %d %d';
        cut_info = textscan(fid, form, 1);
        C = cut_info{4};
    end
end
fclose(fid);

switch abs(ICOMP)
    case 1
        polType = 'linear';
        coorSys = 'spherical';
        E1ff = E1;
        E2ff = E2;
    case 2
        polType = 'circular';
        coorSys = 'spherical';
        E1ff = E2;
        E2ff = E1;
    case 3
        polType = 'linear';
        coorSys = 'Ludwig3';
        E1ff = E2;
        E2ff = E1;
    otherwise
        error(['ICOMP ',num2str(ICOMP),' case not implemented yet'])
end
gridType = 'PhTh';

%% Build the FF obj
x = deg2rad(ph_deg);
y = deg2rad(th_deg);
E3ff = zeros(size(E1ff));
freqUnit = 'Hz';
Prad = ones(1,nr_freq).*4*pi;
[radEff,freq] = deal(ones(1,nr_freq));

FF = FarField(x,y,E1ff,E2ff,E3ff,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);

end
