classdef reflector
    properties
        surface
        rim
        coor
    end
    
    methods
        function obj = reflector(surface,rim,coor)
            if nargin == 0
                obj.surface = paraboloid();
                obj.rim = ellipticalRim();
                obj.coor = CoordinateSystem();
            else
                obj.surface = surface;
                obj.rim = rim;
                obj.coor = coor;
            end
        end
        
        function [surfPoints,rimPoints,surfPointsMesh] = getPointCloud(obj,Npoints,type)
            % Returns a point cloud of the surface and rim points.  The
            % number of points on the surface is of the order N (not exactly - depends on the grid type), and on
            % the rim = sqrt(N) (default N = 1000).
            % type specifies which type of grid to use: ('cart') | 'polar' |
            % 'polarThin' | 'x0' | 'y0' (x0 and y0 returns a linear grid on x=0 or y=0)
            % surfPointsMesh is used for internal plotting, and contains
            % NaNs in the Z data, and is always of the type 'cart' and in a meshgrid style.
            % The actual point cloud removes these NaN points
            
            if nargin == 1
                Npoints = 5000;
                type = 'cart';
            elseif nargin == 2
                type = 'cart';
            elseif nargin == 3
                if isempty(type)
                    type = 'cart';
                end
            end
            N = round(sqrt(Npoints));
            N = N + mod(N,2) + 1;   % Uneven helps in symmetrical systems
            
            % Global coordinate system needed for base changes
            GC = CoordinateSystem;
            
            if nargout == 3
                % Get the surface points on a mesh
                % Use about four times the requested points to make the plot pretty...
                [X,Y] = obj.generateSurfaceGrid(2*N);
                M = obj.rim.isInRim(X,Y);
                Z = obj.surface.getZ(X,Y);
                Pnt = Pnt3D(X,Y,Z);
                Pnt = Pnt.changeBase(GC,obj.coor);
                surfPointsMesh(:,:,1) = Pnt.x;
                surfPointsMesh(:,:,2) = Pnt.y;
                Zrot = Pnt.z;
                Zrot(M == 0) = NaN;
                surfPointsMesh(:,:,3) = Zrot;
                
            end
            
            switch type
                case 'cart'
                    [X,Y] = obj.generateSurfaceGrid(N);
                    x = X(:).';
                    y = Y(:).';
                case {'polar','polarThin'}
                    thinGrid = numel(type) > 5;
                    % Make a polar grid object
                    % Estimate the power from the number of points
                    expN = round(log2(sqrt((Npoints + thinGrid*Npoints)*5))); % Details estimated from the equations in polarGrid. About half of the points are lost when thinning...
                    PG = polarGrid(obj.rim.halfAxis,expN,thinGrid);
                    x = PG.x + obj.rim.centre(1);
                    y = PG.y + obj.rim.centre(2);
                case 'x0'
                    V = obj.rim.cartRim(Npoints);
                    x = linspace(min(V.x),max(V.x),Npoints);
                    y = zeros(size(x));
                case 'y0'
                    V = obj.rim.cartRim(Npoints);
                    y = linspace(min(V.y),max(V.y),Npoints);
                    x = zeros(size(y));
                otherwise
                    error(['Unknown type: ', type])
            end
            M = obj.rim.isInRim(x,y);
            z = obj.surface.getZ(x,y);
            % Remove the points outside the rim
            x(M == 0) = [];
            y(M == 0) = [];
            z(M == 0) = [];
            surfPoints = Pnt3D(x,y,z);
            surfPoints = surfPoints.changeBase(GC,obj.coor);
            
            % Get the actual rim
            Rxy = obj.rim.cartRim(N);
            Rz = obj.surface.getZ(Rxy.x,Rxy.y);
            rimPoints = Pnt3D(Rxy.x,Rxy.y,Rz);
            rimPoints = rimPoints.changeBase(GC,obj.coor);
        end
        
        %% Masking and ray tracing
        function [M,ph_in,th_in] = getMask(obj,coorIn,A,debugPlot)
            % Returns a logical vector M, indicating if the angles in A,
            % defined in the coordinate system coorIn, is pointing towards
            % the reflector or not.
            % A can be a matrix of [ph,th] pairs, or a FarField object.  If
            % it is a FarField object it will be converted to a PhTh grid,
            % and those angles will be used
            
            if nargin < 4
                debugPlot = 0;
            end
            
            if ~isa(A,'FarField')
                [Nr,Nc] = size(A);
                if Nr == 2 && Nc ~= 2
                    A = A.';
                elseif Nc ~= 2
                    error('Unknown input format for A - must be a 2 column matrix or a FarField object');
                end
                FF = FarField(A(:,1),A(:,2),ones(size(A,1),1),zeros(size(A,1),1));
            else
                FF = A;
            end
            ph_in = FF.ph;
            th_in = FF.th;
            
            ph_in = ph_in(:).';
            th_in = th_in(:).';
            
            Npoints = 500^2 + 1;   % Use plenty of points to resolve cases where we're pointing very close to the edge...
            % Get the rim points in 3D space
            [~,rimPoints] = obj.getPointCloud(Npoints);
            rimPointsInCoorIn = changeBase(rimPoints,coorIn);
            ph = rimPointsInCoorIn.ph;
            th = rimPointsInCoorIn.th;
            
            % Interpolate around poles if they are found...
            %             phD = diff(ph);
            %             iJump = find(abs(phD) > 0.99*pi);
            %             for jj = 1:length(iJump)
            %                 if abs(th(iJump) - pi) < deg2rad(0.05) || abs(th(iJump) - 0) < deg2rad(0.05)
            %                     phIn = interp1([iJump:iJump+1],ph(iJump:iJump+1),linspace(iJump,iJump+1));
            %                     phNew = [ph(1:iJump),phIn,ph(iJump+1:end)];
            %                     thIn = interp1([iJump:iJump+1],th(iJump:iJump+1),linspace(iJump,iJump+1));
            %                     thNew = [th(1:iJump),thIn,th(iJump+1:end)];
            %                     ph = phNew;
            %                     th = thNew;
            %                     phD = diff(ph);
            %                     iJump = find(abs(phD) > 0.99*pi);
            %                 end
            %             end
            
            
            iJump = find(abs(th - pi) < deg2rad(0.01));
            if numel(iJump) > 0
                phIn = interp1([iJump-1,iJump+1],[ph(iJump-1),ph(iJump+1)],linspace(iJump-1,iJump+1));
                phNew = [ph(1:iJump-1),phIn,ph(iJump+1:end)];
                thIn = interp1([iJump-1,iJump+1],[th(iJump-1),th(iJump+1)],linspace(iJump-1,iJump+1));
                thNew = [th(1:iJump-1),thIn,th(iJump+1:end)];
                ph = phNew;
                th = thNew;
            end
            
            
            % Test in the TrueView plane - should always make a closed
            % polygon unless at really specific pointings (orthogonal to rim)
            [u,v,w] = PhTh2DirCos(ph,th);
            [Xg,Yg] = DirCos2TrueView(u,v,w);
            [uIn,vIn,wIn] = PhTh2DirCos(ph_in,th_in);
            [XgIn,YgIn] = DirCos2TrueView(uIn,vIn,wIn);
            M = inpolygon(XgIn,YgIn,Xg,Yg);
            % Do final check to determine if we are in fact looking away
            % from the surface...
            % Looking into, or away from, rim.  If [0,0] is in the rim, the
            % curve will close on itself around the origin.  When pointing
            % away, the mask is defined by the region outside the polygon.
            if inpolygon(0,0,Xg,Yg)
                % Estimate the midpoint of the reflector
                RimC = obj.rim.centre;
                P0z = obj.surface.getZ(RimC(1),RimC(2));
                P0 = Pnt3D(RimC(1),RimC(2),P0z);
                % In the global coordinate system
                P0 = P0.changeBase(CoordinateSystem,obj.coor);
                % Get the angle between the input coor z_axis and this
                % point
                V0 = P0-coorIn.origin;
                V0 = [V0.x,V0.y,V0.z]./V0.r;
                psiForward = acos(dot(V0.',coorIn.z_axis));
                psiBackward = acos(dot(V0.',-coorIn.z_axis));
                if psiBackward < psiForward
                    M = ~M;
                end
            end
            
            if debugPlot
                figure
                plot(rad2deg(ph),rad2deg(th),'.-k'), grid on, hold on
                plot(rad2deg(ph_in(M)),rad2deg(th_in(M)),'b.')
                plot(rad2deg(ph_in(~M)),rad2deg(th_in(~M)),'r.')
                figure
                plot(Xg,Yg,'.-k'), grid on, hold on
                plot(XgIn(M),YgIn(M),'b.')
                plot(XgIn(~M),YgIn(~M),'r.')
                axis([-pi,pi,-pi,pi])
            end
        end
        
        function [Pintercept,M] = getRayInterceptPoint(obj,coorIn,ph_in,th_in,NpointsTest,debugPlot)
            % Return the point on the reflector where the ray pointing in the direction
            % [ph_in,th_in] from the coorIn coordinate system will intercept.
            % ph_in and th_in can be vectors of equal length
            % Pintercept is a Pnt3D object in the global coordinate system.
            % M is the associated input ray mask
            % NpointsTest is an optional parameter which determines the
            % resolution on the reflector to test with - keep under 1000 to
            % use to more accurate 3D scettered interpolation scheme
            
            % Just work with internal points, not the rim.  We also only
            % caulculate the nearest internal points for actual valid rays
            % that intercepts the reflector - not those missing it
            
            assert(all(size(ph_in)==size(th_in)),'ph_in and th_in should be the same size');
            if nargin < 5
                NpointsTest = 500;
                debugPlot = 0;
            elseif nargin == 5
                debugPlot = 0;
            end
            ph_in = ph_in(:).';
            th_in = th_in(:).';
            % Get the valid output rays
            M = obj.getMask(coorIn,[ph_in(:),th_in(:)]);
            assert(~any(isnan(M)),'No ray interception points found - coordinate system not pointing at dish')
            ph_out = ph_in(M);
            th_out = th_in(M);
            Nvalid = sum(M);
            
            % Project the reflector surface and rim onto a sphere around the
            % coorIn coordinate system
            [intPoints] = obj.getPointCloud(NpointsTest,'cart');
            % Get points (currently in Global coordinate base) in the coordinate System base
            intPointsInCoorIn = intPoints.changeBase(coorIn);
            % Project onto unit sphere
            int1Sph = Pnt3D.sph(intPointsInCoorIn.ph,intPointsInCoorIn.th,1);
            % Also get the pointsof interest on the unit sphere
            [xInSph,yInSph,zInSph] = sph2cart(ph_in,pi/2-th_in,1);
            in1Sph = Pnt3D(xInSph,yInSph,zInSph);
            
            if NpointsTest <= 1001
                % Use a 3D interpolant on the unit sphere to get the radius
                % to the interception point - can be very expensive for a
                % densely sampled reflector
                R = scatteredInterpolant([int1Sph.x;int1Sph.y;int1Sph.z].',intPointsInCoorIn.r.','natural');
                rInterpolate = R([in1Sph.x(M);in1Sph.y(M);in1Sph.z(M)].').';
                Pintercept = Pnt3D.sph(ph_out,th_out,rInterpolate);
                Pintercept = Pintercept.changeBase(CoordinateSystem(),coorIn);
            else
                % Just get the closest 5 points on the unit sphere and
                % average the correspong results on the reflector
                % Calculate the distances from each input point to each internal
                % and each rim point
                Nint = length(int1Sph.x);
                [xInt,yInt,zInt] = deal(zeros(Nint,Nvalid));
                valIndex = find(M);
                for ii = 1:sum(M)
                    Dint = distanceCart(in1Sph.getNpts(valIndex(ii)),int1Sph).';
                    [~,IintSort] = sort(Dint.');
                    xInt(:,ii) = intPoints.x(IintSort);
                    yInt(:,ii) = intPoints.y(IintSort);
                    zInt(:,ii) = intPoints.z(IintSort);
                end
                %                 pIntClosest = Pnt3D(xInt(1,:),yInt(1,:),zInt(1,:));
                nMean = 5; % Number of closest points to consider
                Pintercept = Pnt3D(mean(xInt(1:nMean,:)),mean(yInt(1:nMean,:)),mean(zInt(1:nMean,:)));
            end
            
            if debugPlot
                int1Sph.plot('marker','.','markerEdgeColor','b','markerSize',5)
                hold on
                in1Sph.plot('marker','.','markerEdgeColor','k','markerSize',5)
                
                figure
                obj.plotMask(coorIn,[ph_out,th_out],7)
                intPoints.plot('marker','.','markerEdgeColor','b','markerSize',5)
                hold on
                %                 pIntClosest.plot('marker','*','markerEdgeColor','g','markerSize',5)
                Pintercept.plot('marker','*','markerEdgeColor','r','markerSize',5)
            end
        end
        
        function [interceptPnt,reflectDir,M] = reflectRays(obj,coorIn,ph_in,th_in,NpointsTest,debugPlot)
            % Returns the point on the reflector, as well as direction of
            % the reflected ray, for the ray pointing in the direction
            % [ph_in,th_in].
            % ph_in and th_in can be vectors of equal length
            % reflectPoint is a Pnt3D object in the global coordinate
            % system, and reflectDir is a [3xNrefl] matrix of unit vectors,
            % in the global coordinate system, of the Nrefl reflected rays.
            % M is the associated input ray mask
            % NpointsTest is an optional parameter which determines the
            % resolution on the reflector to test with - keep under 1000 to
            % use to more accurate 3D scettered interpolation scheme
            
            assert(all(size(ph_in)==size(th_in)),'ph_in and th_in should be the same size');
            if nargin < 5
                NpointsTest = 500;
                debugPlot = 0;
            elseif nargin == 5
                debugPlot = 0;
            end
            ph_in = ph_in(:).';
            th_in = th_in(:).';
            
            % Get the interception points in the global coordinates
            [interceptPnt,M] = getRayInterceptPoint(obj,coorIn,ph_in,th_in,NpointsTest);
            % Get the points in the reflector base coordinate system
            interceptPnt_base = interceptPnt.changeBase(obj.coor);
            % Get the normal vectors here
            Vn = obj.surface.getNorm(interceptPnt_base.x,interceptPnt_base.y);
            % Get the coordinate system in the reflector coordinate system
            coorIn_base = coorIn.redefineToOtherBase(obj.coor);
            % Get the vectors from the coordinate to the reflection points
            Vin = interceptPnt_base.pointMatrix - coorIn_base.origin.pointMatrix;
            VinMag = sqrt(sum(Vin.^2));
            VinN = bsxfun(@rdivide,Vin,VinMag);
            % Reflect the incoming ray
            VoutN = VinN - bsxfun(@times,Vn,2*(dot(Vn,VinN)));
            % Get the points at the tips of the reflected ray for the base
            % change
            reflPointMat_base = interceptPnt_base.pointMatrix + VoutN;
            reflPoint_base = Pnt3D(reflPointMat_base(1,:),reflPointMat_base(2,:),reflPointMat_base(3,:));
            % Get back in global coordinates
            reflPoint = reflPoint_base.changeBase(CoordinateSystem,obj.coor);
            reflectDirTmp = reflPoint - interceptPnt;
            reflectDir = reflectDirTmp.pointMatrix;
            
            
            if debugPlot
                % First plot in the reflector base
                figure
                reflectorPnts = obj.getPointCloud;
                reflectorPnts = reflectorPnts.changeBase(obj.coor);
                reflectorPnts.plot('marker','.','markerEdgeColor','k','markerSize',5)
                hold on
                interceptPnt_base.plot('marker','o','markerEdgeColor','r','markerSize',5)
                interceptPnt_base.plotVect(Vn,'lineColor','r')
                coorIn_base.origin.plotVect(Vin,'lineColor','g')
                interceptPnt_base.plotVect(VoutN,'lineColor','b')
                axis equal
                
                % Now in the actual global coordinate system where we want
                % it
                figure
                obj.plot
                coorIn.plot
                interceptPnt.plot('marker','o','markerEdgeColor','r','markerSize',5)
                plotLines(coorIn.origin,interceptPnt,'lineColor','g')
                interceptPnt.plotVect(reflectDir,'lineColor','b')
            end
        end
        
        
        %% Other functionality and output
        function A = totalArea(obj)
            % ToDo
        end
        
        function writeGRASPsfc(obj,pathName,N)
            % ToDo
        end
        
        function writeGRASPrim(obj,pathName,N)
            % ToDo
        end
        
        %% Plotting
        function plot(obj,N,gridType,surfaceColor)
            % If gridType is empty, or not specified, only the surface is
            % plotted. Otherwise the point cloud is also plotted.
            % See getPointCloud for details on allowable types
            
            if nargin == 1
                N = 10000;
                gridType = [];
                surfaceColor = [0.5,0.5,0.5];
            elseif nargin == 2
                gridType = [];
                surfaceColor = [0.5,0.5,0.5];
            elseif nargin == 3
                surfaceColor = [0.5,0.5,0.5];
            end
            edgeColor = 1.2*surfaceColor;
            rimColor = edgeColor;
            rimWidth = 2;
            % Plot the surface
            [surfPoints,rimPoints,surfPointsMesh] = obj.getPointCloud(N,gridType);
            surf(surfPointsMesh(:,:,1),surfPointsMesh(:,:,2),surfPointsMesh(:,:,3),'FaceColor',surfaceColor,'EdgeColor',edgeColor)
            hold on
            if ~isempty(gridType)
                surfPoints.plot('marker','.','markerEdgeColor','k','markerSize',3)
            end
            axis equal
            % Plot the rim
            %            plot3(rimPoints(1,:),rimPoints(2,:),rimPoints(3,:),'color',rimColor,'LineWidth',rimWidth)
            rimPoints.plot('lineStyle','-','lineWidth',rimWidth,'marker','none','lineColor',rimColor)
            view([140,40])
            xlabel('x-axis (m)')
            ylabel('y-axis (m)')
            zlabel('z-axis (m)')
        end
        
        function plotNorms(obj,Nsurf,Nnorm,scale,gridType)
            % If gridType is empty, or not specified, the normals are
            % plotted on a cartesian grid.
            % See getPointCloud for details on allowable types
            if nargin == 1
                Nsurf = 101;
                Nnorm = Nsurf/5;
                scale = 1;
                gridType = [];
            elseif nargin == 2
                Nnorm = Nsurf/5;
                scale = 1;
                gridType = [];
            elseif nargin == 3
                scale = 1;
            end
            % First just plot the surface
            obj.plot(Nsurf);
            % Now add the normal vectors
            surfPoints = obj.getPointCloud(Nnorm,gridType);
            normalVects = obj.surface.getNorm(surfPoints.x,surfPoints.y).*scale;
            plotVect(surfPoints,normalVects,'lineColor','k')
        end
        
        function plotMask(obj,coorIn,A,scaleRays)
            % Plots a depiction of the rays from coorIn and indicates
            % masking by color.  See getMask for input details
            
            if nargin == 3
                scaleRays = 1;
            end
            [M,ph_in,th_in] = obj.getMask(coorIn,A);
            obj.plot;
            hold on
            coorIn.plot;
            if ~isnan(M)
                hitLen = sum(M);
                misLen = length(ph_in)-sum(M);
                Phit = Pnt3D.sph(ph_in(M),th_in(M),scaleRays);
                Pmis = Pnt3D.sph(ph_in(~M),th_in(~M),scaleRays);
                % Get in the global coordinate system
                GC = CoordinateSystem;
                Phit = Phit.changeBase(GC,coorIn);
                Pmis = Pmis.changeBase(GC,coorIn);
                % Move to origin since we only want relative points for
                % quiver plot
                hitVect = Phit.pointMatrix - coorIn.origin.pointMatrix;
                misVect = Pmis.pointMatrix - coorIn.origin.pointMatrix;
                hitO = repmat(coorIn.origin.pointMatrix,1,hitLen);
                misO = repmat(coorIn.origin.pointMatrix,1,misLen);
                quiver3(hitO(1,:),hitO(2,:),hitO(3,:),hitVect(1,:),hitVect(2,:),hitVect(3,:),'b')
                quiver3(misO(1,:),misO(2,:),misO(3,:),misVect(1,:),misVect(2,:),misVect(3,:),'r')
            end
        end
        
    end
    
    methods (Access = private)
        function [X,Y] = generateSurfaceGrid(obj,N)
            V = obj.rim.cartRim(N);
            V = V.pointMatrix;
            [X,Y] = meshgrid(linspace(min(V(1,:)),max(V(1,:)),N),linspace(min(V(2,:)),max(V(2,:)),N));
        end
    end
end