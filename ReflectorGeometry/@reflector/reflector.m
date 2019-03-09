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
               obj.coor = coordinateSystem();
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
           
           if nargout == 3
               % Get the surface points on a mesh
               % Use about four times the requested points to make the plot pretty...
               [X,Y] = obj.generateSurfaceGrid(2*N);
               M = obj.rim.isInRim(X,Y);
               Z = obj.surface.getZ(X,Y);
               Z(M == 0) = NaN;
               [Xp,Yp,Zp] = changeBaseMeshgrid(obj,X,Y,Z);
               surfPointsMesh(:,:,1) = Xp;
               surfPointsMesh(:,:,2) = Yp;
               surfPointsMesh(:,:,3) = Zp;
           end
           
           switch type
               case 'cart'
                   [X,Y] = obj.generateSurfaceGrid(N);
                   x = X(:);
                   y = Y(:);
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
                   x = linspace(min(V(1,:)),max(V(1,:)),Npoints);
                   y = zeros(size(x));
               case 'y0'
                   V = obj.rim.cartRim(Npoints);
                   y = linspace(min(V(2,:)),max(V(2,:)),Npoints);
                   x = zeros(size(y));
               otherwise
                   error(['Unknown type: ', type])
           end
           x = x(:);
           y = y(:);
           M = obj.rim.isInRim(x,y);
           z = obj.surface.getZ(x,y);
           z(M == 0) = NaN;
           GC = coordinateSystem();
           xyz_prime = changeBase([x,y,z].',GC,obj.coor);
           surfPoints = xyz_prime;
           % Remove columns with NaN z values and return as a 3xN^2 matrix
           NaNCols = any(isnan(surfPoints));
           surfPoints = surfPoints(:,~NaNCols);
           
           % Get the actual rim
           Rxy = obj.rim.cartRim(N);
           Rz = obj.surface.getZ(Rxy(1,:),Rxy(2,:));
           [rimPoints] = changeBase([Rxy;Rz],GC,obj.coor);
       end
       
       function [th,ph,validGraph] = getMaskFunction(obj,coorIn)
           % Returns th, as a function of ph of the apparent rim position 
           % from the coordinate system coorIn.  
           % This information can be used to determine if a certain 
           % direction of radiation, in the coorIn coordinate
           % system, is intercepted by the reflector by checking for <th at
           % ph.
           % validGraph is a logical which indicates if the z-axis actually
           % points at some part of the disch, and therefore presents a
           % valid graph for this fast method of masking...
           
           Npoints = 500^2;   % Use plenty of points to resolve cases where we're pointing very close to the edge...
           % Get the rim points in 3D space
           [~,rimPoints] = obj.getPointCloud(Npoints);
           % Get rim points (currently in Global coordinate base) in the coordinate System base
           rimPointsInCoorIn = changeBase(rimPoints,coorIn);
           % Calculate the vectors between the origin and the rim points
           [ph,el] = cart2sph(rimPointsInCoorIn(1,:),rimPointsInCoorIn(2,:),rimPointsInCoorIn(3,:));
           th = pi/2-el;
           % Check for valid graph (ph must be monotonic)
           validGraph = all(sign(diff(unwrap(ph)))<0) | all(sign(diff(unwrap(ph)))>0);
           if validGraph
               [ph,Isort] = unique(ph);
               th = th(Isort);
           end
       end
       
       function [M,ph_in,th_in] = getMask(obj,coorIn,A)
           % Returns a logical vector M, indicating if the angles in A,
           % defined in the coordinate system coorIn, is pointing towards
           % the reflector or not.
           % A can be a matrix of [ph,th] pairs, or a FarField object.  If
           % it is a FarField object it will be converted to a PhTh grid, 
           % and those angles will be used 
           
           % Interpret the input
           if isa(A,'FarField')
               FF = A.grid2PhTh;
               FF = FF.setXrange('sym');
               ph_in = FF.x;
               th_in = FF.y;
           else
               [Nr,Nc] = size(A);
               if Nr == 2 && Nc ~= 2
                   A = A.';
               elseif Nc ~= 2
                   error('Unknown input format for A - must be a 2 column matrix or a FarField object');
               end
               ph_in = A(:,1);
               th_in = A(:,2);
%                ph_in(ph_in > pi) = ph_in(ph_in > pi) - 2*pi;
               [x_in,y_in,z_in] = sph2cart(ph_in,pi/2-th_in,1);
               [ph_in,el_in,~] = cart2sph(x_in,y_in,z_in);
               th_in = pi/2-el_in;
           end
           [th,ph,validGraph] = obj.getMaskFunction(coorIn);
           if validGraph
               th_test = interp1(ph,th,ph_in,'linear','extrap');
               M = th_in <= th_test;
           else
               M = NaN;
           end
       end
       
       function [Pinternal,Prim] = getRayClosestPoint(obj,coorIn,ph_in,th_in,NpointsTest)
           % Return the closest point on the reflector (internal point, and
           % rim) past which the ray pointing in the direction
           % [ph_in,th_in] from the coorIn coordinate system will pass.
           % ph_in and th_in can be vectors of equal length
           % Pinternal and Prim are returned as [3xNumel(ph_in)] matrices,
           % with the rows [x;y;z] in the global coordinate system
           % NpointsTest is an optional parameter which determines the
           % resolution on the reflector to test with
           
           % Much still ToDo...
           
           assert(all(size(ph_in)==size(th_in)),'ph_in and th_in should be the same size');
           if nargin < 5
               NpointsTest = 200^2;
           end
           ph_in = ph_in(:).';
           th_in = th_in(:).';
           Nin = numel(ph_in);
           % Project the reflector surface and rim onto a sphere around the
           % coorIn coordinate system
           [intPoints,rimPoints] = obj.getPointCloud(NpointsTest);
           % Get points (currently in Global coordinate base) in the coordinate System base
           intPointsInCoorIn = changeBase(intPoints,coorIn);
           rimPointsInCoorIn = changeBase(rimPoints,coorIn);
           % Project onto unit sphere
           [phInt,elInt,~] = cart2sph(intPointsInCoorIn(1,:),intPointsInCoorIn(2,:),intPointsInCoorIn(3,:));
           [phRim,elRim,~] = cart2sph(rimPointsInCoorIn(1,:),rimPointsInCoorIn(2,:),rimPointsInCoorIn(3,:));
           [xIntSph,yIntSph,zIntSph] = sph2cart(phInt,elInt,1);
           intSph = [xIntSph;yIntSph;zIntSph];
           [xRimSph,yRimSph,zRimSph] = sph2cart(phRim,elRim,1);
           rimSph = [xRimSph;yRimSph;zRimSph];
           % Also get the pointsof interest on the unit sphere
           [xInSph,yInSph,zInSph] = sph2cart(ph_in,pi/2-th_in,1);
           inSph = [xInSph;yInSph;zInSph];
           % Calculate the distances from each input point to each internal
           % and each rim point
%            % ToDo: Update this to great circle distances instead of linear 
%            Nint = length(xIntSph);
%            distInt = zeros(Nint,Nin);
%            for ii = length(Nin)
%                dist = sqrt(sum((inSph(:,ii) - intSph).^2,1));
%                distInt(:,ii) = dist.';
% %                [distIntSort,IintSort] = sort(distInt);
% %                Pinternal = intPoints(:,IintSort)
% 
%            end
%            keyboard
           
           
           debugPlot = 1;
           if debugPlot
               plot3(xIntSph,yIntSph,zIntSph,'.b'), hold on
               plot3(xRimSph,yRimSph,zRimSph,'or')
               plot3(xInSph,yInSph,zInSph,'.k')
           end
           
           
       end
       
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
       function plot(obj,N,gridType)
           % If gridType is empty, or not specified, only the surface is
           % plotted. Otherwise the point cloud is also plotted.
           % See getPointCloud for details on allowable types
           
           if nargin == 1
               N = 10000;
               gridType = [];
           elseif nargin == 2
               gridType = [];
           end
           surfaceColor = [0.5,0.5,0.5];
           edgeColor = 1.2*surfaceColor;
           rimColor = edgeColor;
           rimWidth = 2;
           % Plot the surface
           [surfPoints,rimPoints,surfPointsMesh] = obj.getPointCloud(N,gridType);
           surf(surfPointsMesh(:,:,1),surfPointsMesh(:,:,2),surfPointsMesh(:,:,3),'FaceColor',surfaceColor,'EdgeColor',edgeColor)
           hold on
           if ~isempty(gridType)
               plot3(surfPoints(1,:),surfPoints(2,:),surfPoints(3,:),'.','color','k')
           end
           axis equal
           % Plot the rim
           plot3(rimPoints(1,:),rimPoints(2,:),rimPoints(3,:),'color',rimColor,'LineWidth',rimWidth)
           view([140,40])
           xlabel('x-axis (m)')
           ylabel('y-axis (m)')
           zlabel('z-axis (m)')
       end
       
       function plotNorms(obj,Nsurf,Nnorm,gridType)
           % If gridType is empty, or not specified, the normals are
           % plotted on a cartesian grid.
           % See getPointCloud for details on allowable types
           if nargin == 1
               Nsurf = 101;
               Nnorm = Nsurf/5;
               gridType = [];
           elseif nargin == 2 
               Nnorm = Nsurf/5;
               gridType = [];
           end
           % First just plot the surface
           obj.plot(Nsurf);
           % Now add the normal vectors
           [surfPoints] = obj.getPointCloud(Nnorm,gridType);
           normalVects = obj.surface.getNorm(surfPoints(1,:),surfPoints(2,:));
           x0 = surfPoints(1,:);
           y0 = surfPoints(2,:);
           z0 = surfPoints(3,:);
           x1 = surfPoints(1,:) + normalVects(1,:);
           y1 = surfPoints(2,:) + normalVects(2,:);
           z1 = surfPoints(3,:) + normalVects(3,:);
           plot3([x0;x1],[y0;y1],[z0;z1],'k')
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
               hitRays = zeros(3,sum(M));
               misRays = zeros(3,misLen);
               el = pi/2 - th_in;
               [hitX,hitY,hitZ] = sph2cart(ph_in(M),el(M),scaleRays);
               [misX,misY,misZ] = sph2cart(ph_in(~M),el(~M),scaleRays);
               hitVect = [hitX(:).';hitY(:).';hitZ(:).'];
               misVect = [misX(:).';misY(:).';misZ(:).'];
               % Get in the global coordinate system
               GC = coordinateSystem;
               hitVect = changeBase(hitVect,GC,coorIn) - coorIn.origin;
               misVect = changeBase(misVect,GC,coorIn) - coorIn.origin;
               C = coorIn;
               hitO = repmat(C.origin,1,hitLen);
               misO = repmat(C.origin,1,misLen);
               quiver3(hitO(1,:),hitO(2,:),hitO(3,:),hitVect(1,:),hitVect(2,:),hitVect(3,:),'b')
               quiver3(misO(1,:),misO(2,:),misO(3,:),misVect(1,:),misVect(2,:),misVect(3,:),'r')
           end
       end
       
   end
   
   methods (Access = private)
       function [X,Y] = generateSurfaceGrid(obj,N)
           V = obj.rim.cartRim(N);
           [X,Y] = meshgrid(linspace(min(V(1,:)),max(V(1,:)),N),linspace(min(V(2,:)),max(V(2,:)),N));
       end
       
       function [Xp,Yp,Zp] = changeBaseMeshgrid(obj,X,Y,Z)
           % Returns the meshgrid style variable X,Y,Z, defined in the
           % global coordinate system, as Xp,Yp,Zp defined in the
           % reflector coordinate system
            [Nx,Ny] = size(X);
            xyz = [X(:).';Y(:).';Z(:).'];
            GC = coordinateSystem;
            xyz_prime = changeBase(xyz,GC,obj.coor);
            Xp = reshape(xyz_prime(1,:),Nx,Ny);
            Yp = reshape(xyz_prime(2,:),Nx,Ny);
            Zp = reshape(xyz_prime(3,:),Nx,Ny);
       end
   end
end