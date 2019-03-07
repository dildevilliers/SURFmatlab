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
                   x = PG.x;
                   y = PG.y;
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
           xyz_prime = changeBase([x,y,z].',obj.coor);
           surfPoints = xyz_prime;
           % Remove columns with NaN z values and return as a 3xN^2 matrix
           NaNCols = any(isnan(surfPoints));
           surfPoints = surfPoints(:,~NaNCols);
           
           % Get the actual rim
           Rxy = obj.rim.cartRim;
           Rz = obj.surface.getZ(Rxy(1,:),Rxy(2,:));
           [rimPoints] = changeBase([Rxy;Rz],obj.coor);
       end
       
       function [th,ph,coorCentre] = getMaskFunction(obj,x)
           % Returns a th, as a function of ph, in the coordinate system
           % coorCentre, which points at the apparent centre of the
           % reflector with (with the its z-axis), and is positioned at the
           % point x.  This information can be used to determine easily if
           % a certain direction of radiation, in the coorCentre coordinate
           % system, is intercepted by the reflector by checking for <th at
           % ph.

           % ToDo
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
            xyz_prime = changeBase(xyz,obj.coor);
            Xp = reshape(xyz_prime(1,:),Nx,Ny);
            Yp = reshape(xyz_prime(2,:),Nx,Ny);
            Zp = reshape(xyz_prime(3,:),Nx,Ny);
       end
   end
end