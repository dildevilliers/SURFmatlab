classdef dualReflector
    % Exactly the same as the GRASP dual reflector wizard
    % Use the variable names of Granet 2002 paper
    % Default is the MeerKAT dish
    
    properties (SetAccess = private)
        % Granet 2002 paper for description
        % PR is used instead of MR
        
        Dm = 13.5 
        Lm = 10.0254  
        Ls = 2.419  
        th_e = deg2rad(48.89) 
        th_0 = deg2rad(-63.20)
        beta = deg2rad(45.47) 
        sigma = 1 
        th_ext = deg2rad(50) % Extension angle in (rad)
        symFact_ext = 1 % Symmetry factor of the SR extension. 0 is symmetric, 1 is bottom, and -1 is top extension
        Df = 0.2 % Feed aperture diameter
        
        type % String describing the system
        F
        FoD
        a
        f
        e
        alpha
        h
        th_U
        th_L
        Dsx
        Dsy
        dSRPR
        dFPR
        Ht
        Lt
        C_SR
        Feq     % Equivalent focal length
        % Extreme points on the reflectors and aperture
        P0
        P1
        P1e % P1 after SR extension
        P2
        P2e % P2 after SR extension
        Q0
        Q1
        Q2
        R0
        R1
        R2
        % Physical information
        PR_chordX
        PR_chordY
        SR_chordX
        SR_chordY
        apArea        % Projected aperture area
        % Actual objects
        PR  % Primary reflector
        SR  % Secondary reflector
        feedCoor
        apCoor
    end
    
    methods
        % Constructor is over specified - so problems can occur when called
        % directly.  This is to maintain flexibility with defining
        % symmetric and offset systems in one class.  Use one of the
        % meta-constructors below, which are specified according to the
        % Granet paper options and calls this constructor correctly
        function obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df)
            if nargin >= 3 % Set 6, symmetric, min block
                obj.Dm = Dm;
                obj.Lm = Lm;
                obj.th_e = th_e;
            end
            if nargin >= 4 % Set 1, symmetric, no block condition
                obj.Ls = Ls;
            end
            if nargin >= 5 % Set 8 offset
                obj.th_0 = th_0;
                obj.beta = beta;
            end
            if nargin >= 7
                obj.sigma = sigma;
            end
            if nargin >= 8, obj.th_ext = th_ext; end
            if nargin == 9, obj.symFact_ext = symFact_ext; end
            if nargin >= 10
                obj.Df = Df; 
            else
                obj.Df = obj.Dm/25; % Handle default here
            end
            
            if obj.sigma == 1
                obj.type = ['Gregorian'];
                SRhandle = str2func('ellipsoid');
            elseif obj.sigma == -1
                obj.type = ['Cassegrain'];
                SRhandle = str2func('hyperboloid');
            end

            % Build the geometry
            if obj.th_0 == 0 && obj.beta == 0
                % Symmetrical case
                obj.type = ['Symmetric ',obj.type];
                obj.f = obj.Ls*(obj.sigma.*obj.Dm - 4.*obj.Lm.*tan(obj.th_e/2))./(2.*obj.sigma.*obj.Dm + 8.*obj.Ls.*tan(obj.th_e/2));   % (4)
                obj.F = obj.Lm + 2.*obj.f;  % (5)
                obj.a = obj.Ls - obj.f; % (6)
                obj.Dsx = 4.*(obj.Ls - obj.f)./(1./sin(obj.th_e) + obj.sigma.*(16.*obj.F.^2 + obj.Dm.^2)./(8.*obj.F.*obj.Dm));   % (7)
                obj.e = obj.f./obj.a;
                obj.alpha = 0;
                obj.th_U = -2.*atan(obj.Dm/(4.*obj.F));
            else % Offset Case
                obj.type = ['Offset ',obj.type];
                obj.e = (1 - obj.sigma*sqrt(tan(obj.beta/2)/tan((obj.beta - obj.th_0)/2)))./(1 + obj.sigma*sqrt(tan(obj.beta/2)/tan((obj.beta - obj.th_0)/2)));  % (5)
                obj.alpha = 2*atan((obj.e+1)/(obj.e-1)*tan(obj.beta/2));   % (6)
                obj.a = obj.Ls/(2 + (obj.e^2-1)/(obj.e*cos(obj.beta-obj.th_0) + 1));  % (22)
                obj.f = obj.a*obj.e;    % (15)
                obj.th_U = 2*atan((1+obj.e)/(1-obj.e)*tan((obj.alpha-obj.sigma*obj.th_e)/2)) + obj.beta;   % (3)
                obj.F = obj.Dm/(4*(tan(-obj.th_U/2) - tan(-obj.th_0/2)));   % (26)
            end
            obj.h = 2*obj.F*tan(-obj.th_0/2);   % (23)
            obj.th_L = -2*atan((2*obj.h-obj.Dm)/(4*obj.F)); % (4)
            obj.Dsx = -obj.sigma*obj.a*(((obj.e^2-1)*sin(obj.beta-obj.th_U))/(obj.e*cos(obj.beta-obj.th_U)+1) - ((obj.e^2-1)*sin(obj.beta-obj.th_L))/(obj.e*cos(obj.beta-obj.th_L)+1));    % (27)
            obj.dSRPR = obj.h - obj.Dm/2 - obj.a*(obj.sigma-1)/2*(obj.e^2-1)*sin(obj.th_U)/(obj.e*cos(obj.beta-obj.th_U)+1) + obj.a*(obj.sigma+1)/2*(obj.e^2-1)*sin(obj.th_L)/(obj.e*cos(obj.beta-obj.th_L)+1); % (11)
            obj.dFPR = obj.h - obj.Dm/2 + 2*obj.f*sin(obj.beta);   % (10)
            obj.Lt = -obj.a*((obj.sigma + 1)/2)*(((obj.e^2-1)*cos(obj.th_L))/(obj.e*cos(obj.beta-obj.th_L) + 1)) + obj.a*((obj.sigma - 1)/2)*(((obj.e^2-1)*cos(obj.th_U))/(obj.e*cos(obj.beta-obj.th_U) + 1)) - (2*obj.h - obj.Dm)^2/(16*obj.F) + obj.F;  % (12)
            obj.Ht = obj.h + obj.Dm/2 - obj.a*((obj.sigma - 1)/2)*(((obj.e^2-1)*sin(obj.th_L))/(obj.e*cos(obj.beta-obj.th_L) + 1)) + obj.a*((obj.sigma + 1)/2)*(((obj.e^2-1)*sin(obj.th_U))/(obj.e*cos(obj.beta-obj.th_U) + 1));  % (13)
            if obj.th_0 == 0 && obj.beta == 0
                obj.Ht = obj.Dm;    % Special case for symmetric - Ht gets new meaning
            end
            obj.FoD = cot(obj.th_e/2)/4;
            obj.Feq = obj.FoD.*obj.Dm;
            
            % Feed coordinate system needed later and in global coordinates
            feedCoor = coordinateSystem();
            feedCoor.origin = pnt3D(-2*obj.f*sin(obj.beta),0,-2*obj.f*cos(obj.beta)); 
            obj.feedCoor = feedCoor.rotGRASP([obj.beta+obj.alpha,0,0]);

            % The important points
            obj.Q0 = pnt3D(obj.h,0,obj.h^2/(4*obj.F) - obj.F);
            obj.Q1 = pnt3D(obj.h-obj.Dm/2,0,(2*obj.h - obj.Dm)^2/(16*obj.F) - obj.F);
            obj.Q2 = pnt3D(obj.h+obj.Dm/2,0,(2*obj.h + obj.Dm)^2/(16*obj.F) - obj.F);
            
            Rz = obj.Dm/2;
            obj.R0 = obj.Q0;
            obj.R1 = obj.Q1;
            obj.R2 = obj.Q2;
            obj.R0 = obj.R0.setZ(Rz);
            obj.R1 = obj.R1.setZ(Rz);
            obj.R2 = obj.R2.setZ(Rz);
            
            OP0 = obj.sigma*(2*obj.a-obj.Ls);   % (32)
            OP1 = -obj.sigma*obj.a*(obj.e^2-1)/(obj.e*cos(obj.th_L - obj.beta) + 1); % (33)
            OP2 = -obj.sigma*obj.a*(obj.e^2-1)/(obj.e*cos(obj.th_U - obj.beta) + 1); % (34)
            
            obj.P0 = pnt3D(obj.sigma*OP0*sin(obj.th_0),0,obj.sigma*OP0*cos(obj.th_0));
            obj.P1 = pnt3D(obj.sigma*OP1*sin(obj.th_L),0,obj.sigma*OP1*cos(obj.th_L));
            obj.P2 = pnt3D(obj.sigma*OP2*sin(obj.th_U),0,obj.sigma*OP2*cos(obj.th_U));
            
            % Handle the extension of the SR
            % Distribute the amount of extension angle according to the
            % symmetry factor
            % Limit to 1
            if abs(obj.symFact_ext) > 1
                obj.symFact_ext = sign(obj.symFact_ext);
            end
            th_ext1 = (obj.th_ext.*obj.symFact_ext + 1)/2;
            th_ext2 = (-obj.th_ext.*obj.symFact_ext + 1)/2;
            extRotAng = (th_ext2 - th_ext1)/2;
            al_cone = obj.alpha + extRotAng;
            th_cone = obj.th_e + obj.th_ext/2;
            
            if obj.th_ext ~= 0
                obj.type = ['Extended SR ', obj.type];
                % Use equations on p116 of Granet - get the rim positions
                % in the Granet SR coordinate system
                if obj.sigma == 1
                    ph = deg2rad([0,180]);
                elseif obj.sigma == -1
                    ph = deg2rad([180,0]);
                end
                rho_sr = obj.a*(obj.e^2 - 1)./(obj.e*(-sin(al_cone).*sin(th_cone).*cos(ph) + cos(al_cone).*cos(th_cone)) - 1);   % P117 Granet top right
                Primx = (cos(al_cone).*sin(th_cone).*cos(ph) + sin(al_cone).*cos(th_cone)).*rho_sr;
                Primy = sin(th_cone).*sin(ph).*rho_sr;
                Primz = (-sin(al_cone).*sin(th_cone).*cos(ph) + cos(al_cone).*cos(th_cone)).*rho_sr - 2.*obj.f;
                
                % Extension points P1e and P2e in the Granet SR coordinate system
                PeSR = pnt3D(Primx,Primy,Primz);
                Dsxe = abs(diff(PeSR.x));
                % Get in Global coor
                gC = coordinateSystem;
                sC = gC.rotGRASP([obj.beta,0,0]);
                Pe = PeSR.changeBase(gC,sC);
                obj.P1e = pnt3D(Pe.x(1),Pe.y(1),Pe.z(1));
                obj.P2e = pnt3D(Pe.x(2),Pe.y(2),Pe.z(2));
            else
                obj.P1e = obj.P1;
                obj.P2e = obj.P2;
                Dsxe = obj.Dsx;
            end
            
            % Use the extension cone for this calculation
            Cx = (distanceCart(obj.feedCoor.origin,obj.P1e)*sin(al_cone + obj.sigma.*th_cone) + distanceCart(obj.feedCoor.origin,obj.P2e).*sin(al_cone - obj.sigma.*th_cone))/2; % (38)
            Cy = 0; % (38)
            Cz = obj.a*sqrt(1 + Cx^2/(obj.f^2 - obj.a^2)) - obj.f; % (38)
            obj.C_SR = pnt3D(Cx,Cy,Cz);     % In the Granet SR coordinate system
            
            ph = linspace(0,2*pi,1001);
            obj.Dsy = max((2*obj.a*(obj.e^2-1)*sin(obj.th_e).*sin(ph))./(obj.e.*(-sin(obj.alpha).*sin(obj.th_e).*cos(ph) + cos(obj.alpha).*cos(obj.th_e)) - 1));   % (39)
            Dsye = max((2*obj.a*(obj.e^2-1)*sin(th_cone).*sin(ph))./(obj.e.*(-sin(al_cone).*sin(th_cone).*cos(ph) + cos(al_cone).*cos(th_cone)) - 1));   % (39)

            % Derived geometry
            obj.PR_chordX = distanceCart(obj.Q2,obj.Q1);
            obj.PR_chordY = obj.Dm;
            obj.SR_chordX = distanceCart(obj.P2e,obj.P1e);
            obj.SR_chordY = Dsye;
            obj.apArea = pi*(obj.Dm/2)^2;

            % Build the reflectors and coordinate systems
            obj.SR = reflector;
            if obj.sigma == -1
                % Use the standard stuff in Granet for a Cassegrain
                obj.SR.surface = SRhandle(2*obj.a,2*obj.f);
                obj.SR.rim = ellipticalRim([obj.C_SR.x,obj.C_SR.y],[Dsxe,Dsye]./2);
                SRcoor = coordinateSystem;
                obj.SR.coor = SRcoor.rotGRASP([obj.beta,0,0]);
            elseif obj.sigma == 1
                % Use a different SR coordinate system for Gregorians to
                % allow deep sub-reflectors
                % The SR coordinate system rotated to be orthogonal to the
                % SR rim plane with no extension
                % If it is done to the extension points you get different
                % reflectors for extended and not extended cases

%                 SRxPlane = obj.P1e - obj.P2e; % In global coordinates (extended edge plane)
                SRxPlane = obj.P1 - obj.P2; % In global coordinates (non-extended edge plane)
                obj.SR.coor = coordinateSystem(pnt3D(0,0,0),SRxPlane.pointMatrix./SRxPlane.r,[0;-1;0]);
                
                % Project the extension edge points onto the SRcoor
                P1sr = obj.P1e.changeBase(obj.SR.coor);
                P2sr = obj.P2e.changeBase(obj.SR.coor);
                % Get the midpoint
                SRchordXline = P1sr - P2sr;
                SRmidPoint = P2sr + SRchordXline.scale(0.5);
                obj.SR.rim = ellipticalRim([SRmidPoint.x,SRmidPoint.y],[SRchordXline.x,obj.SR_chordY]./2);

                obj.SR.surface = SRhandle(2*obj.a,2*obj.f,pi/2 + (obj.beta - SRxPlane.th));

%                 % Work with Jamnejad 1980 variable names for clarity - all
%                 % trailing *J
%                 alphaJ = th_cone;
%                 betaJ = pi - al_cone;
%                 % Eqn 3.8
%                 gammaJ = atan((cos(betaJ) + obj.e*cos(alphaJ))/sin(betaJ));
%                 
%                 obj.SR.surface = SRhandle(2*obj.a,2*obj.f,pi/2 - gammaJ);
%                 obj.SR.coor = SRcoor.rotGRASP([gammaJ + obj.beta + pi/2,0,pi]);
%                 % Calculate the midpoint of the SR rim
                
                
            end
            
            obj.PR = reflector;
            obj.PR.surface = paraboloid(pnt3D(0,0,-obj.F),obj.F);
            obj.PR.rim = ellipticalRim([obj.Q0.x,obj.Q0.y],[obj.Dm,obj.Dm]./2);
            obj.PR.coor = coordinateSystem();
            
            obj.apCoor = coordinateSystem(obj.R0);
        end
        
        function [rho,drho_dth] = getThRhoMapping(obj,th)
            % Returns the th->rho mapping and its derivative 
            rho = 2.*obj.Feq.*tan(th./2);
            drho_dth = obj.Feq.*sec(th./2).^2;
        end
        
        function [pAp,pReflPR,pReflSR] = rayTrace(obj,ph_in,th_in)
            % Returns the aperture and reflector points for a general 
            % element ray trace
            
            % Get SR reflection points
            [pReflSR,reflectDirSR] = obj.SR.reflectRays(obj.feedCoor,ph_in,th_in);
            % Get the PR reflection points
            [xPR,yPR,zPR] = deal(zeros(size(pReflSR)));
            reflectDirPR = zeros(size(reflectDirSR));
            % Initialise aperture plane normal vector;
            n = repmat([0;0;1],1,size(reflectDirSR,2));    % Normal vector on aperture plane
            % Loop through the directions of interest
            for rr = 1:length(pReflSR.x)
                directionPoint = pnt3D(reflectDirSR(1,rr),reflectDirSR(2,rr),reflectDirSR(3,rr));
                reflCoor = coordinateSystem;
                reflCoor = reflCoor.rotGRASP([directionPoint.th,directionPoint.ph,0]);
                reflCoor.origin = pnt3D(pReflSR.x(rr),pReflSR.y(rr),pReflSR.z(rr));
                [pPR,reflectDirPR(:,rr)] = obj.PR.reflectRays(reflCoor,0,0,5000); % Just pointing along the axis for this one-by-one calculation
                [xPR(rr),yPR(rr),zPR(rr)] = deal(pPR.x,pPR.y,pPR.z);
            end
            pReflPR = pnt3D(xPR,yPR,zPR);
            % Distance from intersection points to plane
            d = dot((obj.apCoor.origin.pointMatrix - pReflPR.pointMatrix),n)./(dot(reflectDirPR,n));
            reflectDirScale = bsxfun(@times,reflectDirPR,d);
            pAp = pReflPR.addVect(reflectDirScale);
        end
        
        function pathLengthStruct = getPathLength(obj,ph_in,th_in)
            % Calculate path length from:
            % feed to SR: FS
            % from SR to PR: SP
            % from PR to aperture: PA
            % from feed to aperture: FA
            if nargin == 1
                th_in = linspace(-obj.th_e*0.99,obj.th_e*0.99,21);
                ph_in = zeros(size(th_in));
            end
            ph_in = ph_in(:).';
            th_in = th_in(:).';
            [pAp,pReflPR,pReflSR] = rayTrace(obj,ph_in,th_in);
            delFS = pReflSR - obj.feedCoor.origin;
            delSP = pReflPR - pReflSR;
            delPA = pAp - pReflPR;
            pathLengthStruct.ph = ph_in;
            pathLengthStruct.th = th_in;
            pathLengthStruct.FS = delFS.r;
            pathLengthStruct.SP = delSP.r;
            pathLengthStruct.PA = delPA.r;
            pathLengthStruct.FA = pathLengthStruct.FS + pathLengthStruct.SP + pathLengthStruct.PA;
        end
        
        function [FFM_F,MaskPointing,Mint] = getMask(obj,A,refl)
            % Returns the reflector mask (for the PR) as a FarField
            % object. Also returns a matrix of pointing directions, as 
            % coordinate system objects, of the final rays after reflection 
            % through the whole system. These can be used to assign
            % background temperature whan calculating antenna temperature.
            % refl can be 'PR' or 'SR'
            % Finally, the logical Mask is returned - true if being masked.
            % A can be a matrix of [ph,th] pairs, or a FarField object.  If
            % it is a FarField object it will be converted to a PhTh grid,
            % and those angles will be used
            if nargin == 1
                A = FarField;
                refl = 'PR';
            end
            if isa(A,'FarField')
                freq = A.freq;
            else
                freq = 1;
            end
            
            switch refl
                case 'PR'
                    [Mint,ph_in,th_in] = obj.PR.getMask(coordinateSystem,A,0);
                    P = repmat(double(Mint(:)),1,numel(freq));
                    FFM_F = FarField.farFieldFromPowerPattern(ph_in(:),th_in(:),P,freq);
                    % Build the pointing matrix - already at globalCoor base
                    MaskPointing(size(ph_in)) = coordinateSystem;
                    % Fix the non-masked position pointing angles - in the global
                    % coordinate system
                    maskI = find(~Mint);
                    for m0 = maskI
                        MaskPointing(m0) = MaskPointing(m0).rotGRASP([th_in(m0),ph_in(m0),0]);
                    end
                case 'SR'
                    [M,ph_in,th_in] = obj.SR.getMask(obj.feedCoor,A,0);
                    % Also require an SR with no extension
                    if obj.th_ext > 0
                        SRnoExt = dualReflector(obj.Dm,obj.Lm,obj.th_e,obj.Ls,obj.th_0,obj.beta,obj.sigma,0,0);
                        MnoExt = SRnoExt.SR.getMask(obj.feedCoor,A,0);
                    else 
                        MnoExt = M;
                    end
                    Mint = M + MnoExt;  % 0 where no mask, 1 on the extension, 2 in the actual SR
                    % The FarField object is based on the intersected
                    % extension - all the information is here
                    P = repmat(double(Mint(:)),1,numel(freq));
                    FFM_F = FarField.farFieldFromPowerPattern(ph_in(:),th_in(:),P,freq);
                    
                    % Build the pointing matrix
                    MaskPointing(size(ph_in)) = coordinateSystem;
                    % First those outside the SR mask - centered at feed currently
                    [MaskPointing(Mint<=1).base] = deal(obj.feedCoor);
                    % Those in mask - origin at global base and pointing up so do
                    % nothing...
                    % Fix the non-masked position pointing angles - in the global
                    % coordinate system
                    % First the outside ones - no reflection
                    mask0 = find(Mint == 0);
                    for m0 = mask0
                        MaskPointing(m0) = MaskPointing(m0).rotGRASP([th_in(m0),ph_in(m0),0]);
                        MaskPointing(m0) = MaskPointing(m0).getInGlobal;    % Rotate to global Coor
                        MaskPointing(m0).origin = pnt3D;    % Force to centre of global Coor
                    end
                    % And reflect the rays pointing at the extension
                    mask1 = find(Mint == 1);
                    % Reflection directions already in Global Coordinates
                    [~,reflectDir] = obj.SR.reflectRays(obj.feedCoor,ph_in(mask1),th_in(mask1));
                    ii = 0;
                    for m1 = mask1
                        ii = ii+1;
                        directionPoint = pnt3D(reflectDir(1,ii),reflectDir(2,ii),reflectDir(3,ii));
                        MaskPointing(m1) = MaskPointing(m1).rotGRASP([directionPoint.th,directionPoint.ph,0]);
                        MaskPointing(m1).origin = pnt3D;    % Force to centre of global Coor
                    end
                otherwise
                    error(['Undefined input: ', refl]);
            end
        end
        
        %% Plotting
        function plot(obj,N)
            if nargin == 1
                N = 101;
            end
            lineWidthRefl = 2;
            lineWidthRays = 0.5;
            % Plot the reflector surface in the symmetry plane
            [surfPointsPR] = obj.PR.getPointCloud(N,'x0');
            [surfPointsSR] = obj.SR.getPointCloud(N,'x0');
            plot(surfPointsPR.x,surfPointsPR.z,'k','linewidth',lineWidthRefl), hold on, grid on
            plot(surfPointsSR.x,surfPointsSR.z,'r','linewidth',lineWidthRefl)
            % Plot the foci
            globalCoor = coordinateSystem;
            F0 = obj.SR.surface.F0.changeBase(globalCoor,obj.SR.coor);
            F1 = obj.SR.surface.F1.changeBase(globalCoor,obj.SR.coor);
            plot(F0.x,F0.z,'ko')
            plot(F1.x,F1.z,'ko')
            % Overlay plot the original SR if an extension is present
            if obj.th_ext > 0
                objNoExt = dualReflector(obj.Dm,obj.Lm,obj.th_e,obj.Ls,obj.th_0,obj.beta,obj.sigma,0,0);
                [surfPointsSRNoExt] =  objNoExt.SR.getPointCloud(N,'x0');
                plot(surfPointsSRNoExt.x,surfPointsSRNoExt.z,'k','linewidth',lineWidthRefl)
                F0noExt = objNoExt.SR.surface.F0.changeBase(globalCoor,objNoExt.SR.coor);
                F1noExt = objNoExt.SR.surface.F1.changeBase(globalCoor,objNoExt.SR.coor);
                plot(F0noExt.x,F0noExt.z,'ro')
                plot(F1noExt.x,F1noExt.z,'ro')
            end
            % Plot the feed point
            plot(obj.feedCoor.origin.x,obj.feedCoor.origin.z,'k.','markersize',2)
            % Plot the edge rays
            xR = [obj.R2.x,obj.Q2.x,obj.P2.x,obj.feedCoor.origin.x,obj.P1.x,obj.Q1.x,obj.R1.x];
            zR = [obj.R2.z,obj.Q2.z,obj.P2.z,obj.feedCoor.origin.z,obj.P1.z,obj.Q1.z,obj.R1.z];
            plot(xR,zR,'k','linewidth',lineWidthRays)
            xlabel('x (m)')
            ylabel('z (m)')
            axis equal
        end
        
        function plot3D(obj,N,coorFlag)
            % N is the approximated root of the number of points
            % coorFlag is a vector of length 3 indicating if the coordinate
            % systems should be plotted: [feedCoor,apertureCoor,reflectorCoor]
            % Default is [1,0,0]
            if nargin == 1
                N = 10000;
                coorFlag = [1,0,0];
            elseif nargin  == 2
                coorFlag = [1,0,0];
            end
            obj.SR.plot(N), hold on
            obj.PR.plot(N)
            if coorFlag(1), obj.feedCoor.plot(obj.Dm/10); end
            if coorFlag(2), obj.apCoor.plot(obj.Dm/10); end
            if coorFlag(3)
                obj.PR.coor.plot(obj.Dm/10); 
                obj.SR.coor.plot(obj.Dm/10);
            end
            xlabel('x-axis (m)')
            ylabel('y-axis (m)')
            zlabel('z-axis (m)')
            % Plot the original SR rim if an extension is present
            if obj.th_ext > 0
                objNoExt = dualReflector(obj.Dm,obj.Lm,obj.th_e,obj.Ls,obj.th_0,obj.beta,obj.sigma,0,0);
                objNoExt.SR.plot(10000,[],[1,1,1].*0.2)
            end
        end
        
        function plotRayTrace(obj,Nray,Nrefl)
            if nargin == 1
                Nray = 21;
                Nrefl = 101;
            elseif nargin == 2
                Nrefl = 101;
            end
            th_in = linspace(-obj.th_e*0.99,obj.th_e*0.99,Nray);
            ph_in = zeros(size(th_in));
            obj.plot3D(Nrefl)
            hold on
            [pAp,pReflPR,pReflSR] = rayTrace(obj,ph_in,th_in);
            rayColor = ones(1,3).*0;
            rayWidth = 0.5;
            plotLines(obj.feedCoor.origin,pReflSR,'lineColor',rayColor,'lineWidth',rayWidth)
            plotLines(pReflSR,pReflPR,'lineColor',rayColor,'lineWidth',rayWidth)
            plotLines(pReflPR,pAp,'lineColor',rayColor,'lineWidth',rayWidth)
            view([0,0])
        end
        
    end
    
    % Make a ton of different constructor options...
    methods (Static = true)
        function obj = ODR1(Dm,F,h,Dsx,beta,sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption1(Dm, F, h, Dsx, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR2(Dm, F, h, Ls, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption2(Dm, F, h, Ls, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR3(Dm, F, h, DFPR, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption3(Dm, F, h, DFPR, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR4(Dm, F, h, Lt, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption4(Dm, F, h, Lt, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR5(Dm, F, h, Ht, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption5(Dm, F, h, Ht, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR6(Dm, F, h, dSRPR, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption6(Dm, F, h, dSRPR, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR7(Dm, th_0, dFPR, Ls, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption7(Dm, th_0, dFPR, Ls, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR8(Dm, th_0, th_e, Ls, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [~,~,~,~,~,~,~,~,~,~,Lm] = ODRoption8(Dm, th_0, th_e, Ls, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR9(Dm, th_0, th_e, Dsx, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption9(Dm, th_0, th_e, Dsx, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR10(Dm, th_0, th_e, dSRPR, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption10(Dm, th_0, th_e, dSRPR, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR11(Dm, th_0, th_e, Lt, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption11(Dm, th_0, th_e, Lt, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = ODR12(Dm, th_0, th_e, Ht, beta, sigma,th_ext,symFact_ext,Df)
            if nargin <= 6, th_ext = 0; end
            if nargin <= 7, symFact_ext = 1; end
            if nargin <= 8, Df = 0.2; end
            [Dm, th_0, th_e, Ls, Lm, beta] = ODRoption12(Dm, th_0, th_e, Ht, beta, sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        
        function obj = SDR1(Dm,Lm,Ls,th_e,sigma,th_ext,symFact_ext,Df)
            if nargin <= 5, th_ext = 0; end
            if nargin <= 6, symFact_ext = 0; end
            if nargin <= 7, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRoption1(Dm,Lm,Ls,th_e,sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDR2(Dm,F,Lm,th_e,sigma,th_ext,symFact_ext,Df)
            if nargin <= 5, th_ext = 0; end
            if nargin <= 6, symFact_ext = 0; end
            if nargin <= 7, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRoption2(Dm,F,Lm,th_e,sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDR3(Dm,F,Ls,th_e,sigma,th_ext,symFact_ext,Df)
            if nargin <= 5, th_ext = 0; end
            if nargin <= 6, symFact_ext = 0; end
            if nargin <= 7, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRoption3(Dm,F,Ls,th_e,sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDR4(F,Ds,Ls,th_e,sigma,th_ext,symFact_ext,Df)
            if nargin <= 5, th_ext = 0; end
            if nargin <= 6, symFact_ext = 0; end
            if nargin <= 7, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRoption4(F,Ds,Ls,th_e,sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDR5(Lm,Ds,Ls,th_e,sigma,th_ext,symFact_ext,Df)
            if nargin <= 5, th_ext = 0; end
            if nargin <= 6, symFact_ext = 0; end
            if nargin <= 7, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRoption5(Lm,Ds,Ls,th_e,sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDR6(Dm,F,Ds,th_e,sigma,th_ext,symFact_ext,Df)
            if nargin <= 5, th_ext = 0; end
            if nargin <= 6, symFact_ext = 0; end
            if nargin <= 7, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRoption6(Dm,F,Ds,th_e,sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDR7(Dm,Ds,Ls,th_e,sigma,th_ext,symFact_ext,Df)
            if nargin <= 5, th_ext = 0; end
            if nargin <= 6, symFact_ext = 0; end
            if nargin <= 7, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRoption7(Dm,Ds,Ls,th_e,sigma);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDRmb1(Dm,F,Lm,sigma, Df,th_ext,symFact_ext)
            if nargin <= 4, th_ext = 0; end
            if nargin <= 5, symFact_ext = 0; end
            if nargin <= 6, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRmboption1(Dm,F,Lm,sigma, Df);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDRmb2(Dm,F,th_e, sigma, Df,th_ext,symFact_ext)
            if nargin <= 4, th_ext = 0; end
            if nargin <= 5, symFact_ext = 0; end
            if nargin <= 6, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRmboption2(Dm,F,th_e, sigma, Df);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDRmb3(Dm,F,Ds, sigma, Df,th_ext,symFact_ext)
            if nargin <= 4, th_ext = 0; end
            if nargin <= 5, symFact_ext = 0; end
            if nargin <= 6, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRmboption3(Dm,F,Ds, sigma, Df);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDRmb4(Dm,Lm,Ds, sigma, Df,th_ext,symFact_ext)
            if nargin <= 4, th_ext = 0; end
            if nargin <= 5, symFact_ext = 0; end
            if nargin <= 6, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRmboption4(Dm,Lm,Ds, sigma, Df);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDRmb5(Dm,Ds,th_e, sigma, Df,th_ext,symFact_ext)
            if nargin <= 4, th_ext = 0; end
            if nargin <= 5, symFact_ext = 0; end
            if nargin <= 6, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRmboption5(Dm,Ds,th_e, sigma, Df);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDRmb6(Dm,Lm,th_e, sigma, Df,th_ext,symFact_ext)
            if nargin <= 4, th_ext = 0; end
            if nargin <= 5, symFact_ext = 0; end
            if nargin <= 6, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRmboption6(Dm,Lm,th_e, sigma, Df);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
        function obj = SDRmb7(Dm,Ls,th_e, sigma, Df,th_ext,symFact_ext)
            if nargin <= 4, th_ext = 0; end
            if nargin <= 5, symFact_ext = 0; end
            if nargin <= 6, Df = 0.2; end
            [Dm,~,Lm,~,Ls,~,~,th_e] = SDRmboption7(Dm,Ls,th_e, sigma, Df);
            obj = dualReflector(Dm,Lm,th_e,Ls,th_0,beta,sigma,th_ext,symFact_ext,Df);
        end
    end
    
end







































