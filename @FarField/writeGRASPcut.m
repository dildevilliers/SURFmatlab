function [] = writeGRASPcut(obj,pathName)

% function [] = writeGRASPcut(obj,pathName)
% Function to write GRASP cut files to the location PathName 
% (.cut appended automatically)
% th must be sampled uniformly 

% Dirk de Villiers
% Created: 2019-04-22

assert(obj.isGridUniform,'Standard uniform grid required for GRASP cut write')
assert(strcmp(obj.gridTypeBase,'PhTh'),'Base grid must be PhTh for GRASP cut write')

% Sort out the data
th_vect = unique(obj.th);
ph_vect = unique(obj.ph);
Nth = obj.Ny;
Nph = obj.Nx;

thD = mean(diff(th_vect));

% Set up file variables
HEADER = 'Cut file generated from MATLAB FarField object';
V_INI = rad2deg(th_vect(1));
V_INC = rad2deg(thD);
V_NUM = Nth;
if strcmp(obj.polType,'linear')
    if strcmp(obj.coorType,'spherical')
        ICOMP = 1;  % linear th ph polarization...
        E1real = real(obj.E1);
        E1imag = imag(obj.E1);
        E2real = real(obj.E2);
        E2imag = imag(obj.E2);
    elseif strcmp(obj.coorType,'Ludwig3')
        ICOMP = 3;  % linear Ludwig3 polarization...
        E1real = real(obj.E2);
        E1imag = imag(obj.E2);
        E2real = real(obj.E1);
        E2imag = imag(obj.E1);
    end
elseif strcmp(obj.polType,'circular')
    ICOMP = 2;  % circular polarization...
    E1real = real(obj.E2);
    E1imag = imag(obj.E2);
    E2real = real(obj.E1);
    E2imag = imag(obj.E1);
end
ICUT = 1; % Standard polar cut where phi is fixed (nr_cuts) and th is varied
NCOMP = 2;  % for farfields only 2 components needed

% Create and open the file
PN = [pathName,'.cut'];
fid = fopen(PN,'wt');

% Check for identical start and end ph values
if abs((ph_vect(end)-2*pi) - ph_vect(1)) < eps
    phEnd = Nph-1; 
else
    phEnd = Nph; 
end

for ff = 1:obj.Nf
    for pp = 1:phEnd
        C = rad2deg(ph_vect(pp));
        
        % Header lines
        form = '%s\n';
        fprintf(fid,form,[HEADER]);
        
        form = '%.10E %.10E %i %.10E %i %i %i\n';
        fprintf(fid,form,[V_INI, V_INC, V_NUM, C, ICOMP, ICUT, NCOMP]);
        
        startPos = (pp-1)*Nth + 1;
        stopPos = pp*Nth;
        E1r = E1real(startPos:stopPos,ff);
        E1i = E1imag(startPos:stopPos,ff);
        E2r = E2real(startPos:stopPos,ff);
        E2i = E2imag(startPos:stopPos,ff);
        
        form = '%.10E\t%.10E\t%.10E\t%.10E\n';
        for tt = 1:Nth
            fprintf(fid,form,[E1r(tt), E1i(tt), E2r(tt), E2i(tt)]);
        end
    end
end

fclose(fid);
