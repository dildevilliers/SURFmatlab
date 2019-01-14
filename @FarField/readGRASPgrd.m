function FF = readGRASPgrd(pathName)

% function FF = readGRASPgrd(filePathName)
% Reads a GRASP .grd file and returns the FarField object.
% Not all GRASP functionality supported yet...

[fid] = fopen([pathName,'.grd']);

E1 = [];
E2 = [];
E3 = [];

freqMarker = 'FREQUENCIES';
startMarker = '++++';

% Read the field header
while 1
    a = fgetl(fid);
    if strncmp(a,freqMarker,11) % Read the frequencies info
        freqUnit = regexp(a, '(?<=\[)[^)]*(?=\])', 'match', 'once');
        % Keep reading lines until all frequencies read
        freq = [];
        while 1
            a = fgetl(fid);
            if strcmp(a,startMarker), break; end
            freq = [freq,str2num(a)];
        end
    end
    if strcmp(a,startMarker)
        a = fgetl(fid);
        a = fgetl(fid);
        fieldInfo = str2num(a);
        NSET = fieldInfo(1);
        ICOMP = fieldInfo(2);
        NCOMP = fieldInfo(3);
        IGRID = fieldInfo(4);
        [IX,IY] = deal(zeros(1,NSET));
        for ff = 1:NSET
            a = fgetl(fid);
            centreInfo = str2num(a);
            IX(ff) = centreInfo(1);
            IY(ff) = centreInfo(2);
        end
        break;
    end
end
% Check if all beams have the same grid
cI1 = centreInfo(1,:);
comp = any(bsxfun(@minus,centreInfo,cI1),2);
if sum(comp) > 0
    error('All the field sets must have the same grid - here there are different centre points...');
end

% Front matter done - read the NSET frequency blocks
for ff = 1:NSET
    a = fgetl(fid);
    gridInfo1 = str2num(a);
    XS = gridInfo1(1);
    YS = gridInfo1(2);
    XE = gridInfo1(3);
    YE = gridInfo1(4);
    a = fgetl(fid);
    gridInfo2 = str2num(a);
    NX = gridInfo2(1);
    NY = gridInfo2(2);
    KLIMIT = gridInfo2(3);
    if KLIMIT == 1
        %ToDo
        error('KLIMIT = 1 not implemented yet...')
    else
        IS = 1;
        JS = 1;
        IE = NX;
        JE = NY;
    end
    if ff == 1  % Just build the grid once - assume they are all the same (ToDo: build a check later)
        DX = (XE-XS)/(NX - 1);
        DY = (YE-YS)/(NY - 1);
        XCEN = DX*IX(ff);
        YCEN = DY*IY(ff);
        X = XCEN + XS+DX.*((IS:IE) - 1);
        Y = YCEN + YS+DY.*((JS:JE) - 1);
    end
    
    if NCOMP == 2
        form = '%f %f %f %f';
        fieldData = textscan(fid, form, NX*NY);
        E1 = [E1,(fieldData{1} + 1i.*fieldData{2})];
        E2 = [E2,(fieldData{3} + 1i.*fieldData{4})];
    elseif NCOMP == 3
        form = '%f %f %f %f %f %f';
        fieldData = textscan(fid, form, NX*NY);
        E1 = [E1,(fieldData{1} + 1i.*fieldData{2})];
        E2 = [E2,(fieldData{3} + 1i.*fieldData{4})];
        E3 = [E3,(fieldData{5} + 1i.*fieldData{6})];
    end
    % Dummy read
    a = fgetl(fid);
end
fclose(fid);

%% Build the object
[Xmat,Ymat] = ndgrid(X,Y);
x = Xmat(:);
y = Ymat(:);
switch ICOMP
    case 1
        polType = 'linear';
        coorSys = 'spherical';
    case 2
        polType = 'circular';
        coorSys = 'spherical';
    case 3
        polType = 'linear';
        coorSys = 'Ludwig3';
    otherwise
        error(['ICOMP ',num2str(ICOMP),' case not implemented yet'])
end
switch IGRID
    case 1
        gridType = 'UV';
    case 4
        gridType = 'AzEl';
    case 5
        gridType = 'TrueView';
    case 6 
        gridType = 'ElAz';
    case 7
        gridType = 'PhTh';
    otherwise
        error(['IGRID ',num2str(IGRID),' case not implemented yet'])
end
if isempty(E3)
    E3 = zeros(size(E1));
end
% keyboard;
Prad = ones(size(freq)).*4*pi;
radEff = ones(size(freq));
FF = FarField(x,y,E1,E2,E3,freq,Prad,radEff,coorSys,polType,gridType,freqUnit);
% FF.polType = polType;
% FF.coorSys = coorSys;
% FF.gridType = gridType;
% FF.freqUnit = freqUnit;
FF = setEnames(FF);
FF = setXYnames(FF);
FF = setPhTh(FF);
FF = setFreq(FF);
FF = setBase(FF);


