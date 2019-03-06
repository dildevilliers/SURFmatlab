classdef polarGrid
    % Makes a (thinned) polar grid as described in the 1990 shaping paper
    % by Kildal in TAP
    
    properties
        rhoMax(1,1) double {mustBeReal, mustBeFinite} = 1 % radius in (m)
        expN(1,1) double {mustBeReal, mustBeFinite} = 5   % number of spokes = 2^expN
        thinFlag(1,1) logical = false
    end
    
    properties (SetAccess = private)
        N % Number of outer spokes
        M % Number of rings
        N_m % Number of spokes as a function of M
        delta_n % Thinning factor as function of M
        rho     % vector of rho values
        phi     % vector of phi values
        x       % vector fo x values
        y       % vector of y values
    end
    
    methods
        function obj = polarGrid(rM,I,thinGrid)
            if nargin >= 1
                obj.rhoMax = rM;
                if nargin >= 2
                    obj.expN = I;
                    if nargin == 3
                        obj.thinFlag = thinGrid;
                    end
                end
            end
            
            if obj.expN >= 11
                reply = input(['Are you sure that expN = ',num2str(obj.expN),'? This will make a huge grid! Y/N [N]:'],'s');
                if isempty(reply)
                    reply = 'N';
                end
                if strcmp(reply,'N')
                    return;
                end
            end
            % Build the grid thinning information
            obj.N = 2^obj.expN;
            obj.M = 1 + 3*obj.N/16;
            
            i_m = ones(1,obj.M);
            if obj.thinFlag
                i_m(1) = 3;
            else
                i_m(1) = obj.expN;
            end
            for mm = 2:obj.M
                if obj.thinFlag
                    if mm <= obj.M
                        i_m(mm) = round(log10(16*mm/3)/log10(2));
%                     else % Keep number of phi points the same for the outer th angles
%                         i_m(mm) = i_m(obj.M);
                    end
                else
                    i_m(mm) = obj.expN;
                end
            end
            obj.N_m = 2.^i_m;
            obj.delta_n = obj.N./obj.N_m;
            
            % Build the actual grid
            Npoints = sum(obj.N_m);
            [obj.phi,obj.rho] = deal(zeros(Npoints,1));
            
            rho_grid_delta = obj.rhoMax/(obj.M-1);
            
            obj.rho(1:obj.N_m(1)) = 0*rho_grid_delta;
            obj.phi(1:obj.N_m(1)) = linspace(0,2*pi,obj.N_m(1));
            for mm = 2:obj.M
                obj.rho(sum(obj.N_m(1:mm-1))+1:sum(obj.N_m(1:mm))) = (mm-1)*rho_grid_delta;
                ph_delta = 2*pi/obj.N_m(mm);
                obj.phi(sum(obj.N_m(1:mm-1))+1:sum(obj.N_m(1:mm))) = 0:ph_delta:(obj.N_m(mm)-1)*ph_delta;
            end
            obj.x = obj.rho.*cos(obj.phi);
            obj.y = obj.rho.*sin(obj.phi);
        end
        
        function polar(obj)
            polar(obj.phi,obj.rho,'k.')
        end
        
        function plot(obj)
            plot(obj.x,obj.y,'k.')
            axis equal
            grid on
        end
        
    end
end


