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
       
       function [surfPoints,rimPoints,surfPointsMesh] = getPointCloud(obj,N)
           % Returns a point cloud of the surface and rim points.  The
           % number of points on the surface is of the order N^2, and on
           % the rim = N (default N = 101).
           % The mesh version is used for internal plotting, and contains
           % NaNs in the Z data.  The actual point cloud removes these
           % points
           
           if nargin == 1
               N = 101;
           end
           
           % Get the surface points on a mesh
           [X,Y] = obj.generateSurfaceGrid(N);
           M = obj.rim.isInRim(X,Y);
           Z = obj.surface.getZ(X,Y);
           Z(M == 0) = NaN;
           [Xp,Yp,Zp] = changeBaseMeshgrid(obj,X,Y,Z);
           surfPointsMesh(:,:,1) = Xp;
           surfPointsMesh(:,:,2) = Yp;
           surfPointsMesh(:,:,3) = Zp;
           % Remove columns with NaN z values and return as a 3xN^2 matrix
           surfPoints = [Xp(:).';Yp(:).';Zp(:).'];
           NaNCols = any(isnan(surfPoints));
           surfPoints = surfPoints(:,~NaNCols);
           
           % Get the actual rim
           Rxy = obj.rim.cartRim;
           Rz = obj.surface.getZ(Rxy(1,:),Rxy(2,:));
           [rimPoints] = changeBase([Rxy;Rz],obj.coor);
       end
       
       function B = doesIntersect(obj,coorOfSource)
           % Checks if the z-axis of the coordinate system in coorOfSource
           % intersects the reflector object - just needs to be in the rim
           % in the plane of the rim
       end
       
       function A = projectedArea(obj)
           
       end
       
       function A = totalArea(obj)
           
       end
       
       function writeGRASPsfc(obj,pathName,N)
           
       end
       
       function writeGRASPrim(obj,pathName,N)
           
       end

       %% Plotting
       function plot(obj,N)
           if nargin == 1
               N = 101;
           end
           surfaceColor = [0.5,0.5,0.5];
           edgeColor = 1.2*surfaceColor;
           rimColor = edgeColor;
           rimWidth = 2;
           % Plot the surface
           [surfPoints,rimPoints,surfPointsMesh] = getPointCloud(obj,N);
           surf(surfPointsMesh(:,:,1),surfPointsMesh(:,:,2),surfPointsMesh(:,:,3),'FaceColor',surfaceColor,'EdgeColor',edgeColor)
           hold on
%            plot3(surfPoints(1,:),surfPoints(2,:),surfPoints(3,:),'.','color',edgeColor)
           axis equal
           % Plot the rim
           plot3(rimPoints(1,:),rimPoints(2,:),rimPoints(3,:),'color',rimColor,'LineWidth',rimWidth)
       end
       
   end
   
   methods (Access = private)
       function [X,Y] = generateSurfaceGrid(obj,N)
           V = obj.rim.cartRim(N);
           [X,Y] = meshgrid(linspace(min(V(1,:)),max(V(1,:)),N),linspace(min(V(2,:)),max(V(2,:)),N));
       end
       
       function [Xp,Yp,Zp] = changeBaseMeshgrid(obj,X,Y,Z)
           % Returns the meshgrid style variable X,Y,Z, defined in the
           % global coordinate system, as Xp,Yp, and Zp defined in the
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