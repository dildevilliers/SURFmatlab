function obj1 = fastPhTh2AstroGrid(obj,ngt,xn,yn)

%assuming in this first iteration simply that we want to convert PhTh to
%GalLongLat

assert(strcmp(obj.gridTypeBase,'PhTh'),'FarField base grid type must be PhTh');
obj = obj.grid2AzAlt;

switch ngt
    case 'AzAlt'
        error('PhTh to AzAlt not yet implemented')
    case 'RAdec'
        error('PhTh to RAdec not yet implemented')
    case 'GalLongLat'
        obj1 = obj.grid2AzAlt;
        obj1 = obj1.grid2RAdec;
        obj1 = obj1.grid2GalLongLat;
        U = obj1.getU;
        for ff = 1:obj.Nf
            scatint = scatteredInterpolant([obj1.x obj1.y],U(:,ff));
            Un(:,ff) = scatint(xn,yn);
        end
        obj1 = FarField.farFieldFromPowerPattern(xn,yn,Un,obj.freq,[],[],'GalLongLat');
    otherwise error('newGridType string not recognised')
end



end