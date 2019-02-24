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
       
       function [th,ph,coorCentre] = getMaskFunction(obj,x)
           % Returns a th, as a function of ph, in the coordinate system
           % coorCentre, which points at the apparent centre of the
           % reflector with (with the its z-axis), and is positioned at the
           % point x.  This information can be used to determine easily if
           % a certain direction of radiation, in the coorCentre coordinate
           % system, is intercepted by the reflector by checking for <th at
           % ph.

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
           view([140,40])
           xlabel('x-axis (m)')
           ylabel('y-axis (m)')
           zlabel('z-axis (m)')
       end
       
       function plotNorms(obj,Nsurf,Nnorm)
           if nargin == 1
               Nsurf = 101;
               Nnorm = Nsurf/5;
           elseif nargin == 2 
               Nnorm = Nsurf/5;
           end
           % First just plot the surface
           obj.plot(Nsurf);
           % Now add the normal vectors
           [surfPoints] = getPointCloud(obj,Nnorm);
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