classdef coordinateSystem
    properties
      origin(3,1) double {mustBeReal, mustBeFinite} = [0;0;0] % [x;y;z] in (m) 
      x_axis(3,1) double {mustBeReal, mustBeFinite} = [1;0;0] % x-axis direction
      y_axis(3,1) double {mustBeReal, mustBeFinite} = [0;1;0] % y-axis direction
      base = [] % Can be another coordinate system object, or, if empty, is assumed to be the global coordinate system 
    end
   
    properties (SetAccess = private)
       z_axis 
    end
   methods
       function obj = coordinateSystem(origin,x_axis,y_axis,base)
           if nargin == 0
           elseif nargin == 1
               obj.origin = origin;
           elseif nargin == 2
               obj.origin = origin;
               obj.x_axis = x_axis;
           elseif nargin == 3
               obj.origin = origin;
               obj.x_axis = x_axis;
               obj.y_axis = y_axis;
           else
               obj.origin = origin;
               obj.x_axis = x_axis;
               obj.y_axis = y_axis;
               obj.base = base;
           end
           % Check for valid inputs
           nX = norm(obj.x_axis);
           nY = norm(obj.y_axis);
           z = cross(obj.x_axis,obj.y_axis)/(nX*nY);
           if abs(norm(z)-1) > 1e-10
               error('x_axis and y_axis must be orthogonal');
           end
           % Set the z-axis direction
           obj = obj.setZaxis;
           % Normalise the unit vectors
           obj = obj.normAxis;
       end
       
       function obj = set2Base(obj)
           if class(obj.base) == 'coordinateSystem'
               obj = obj.base;
           end
       end
       
       %% Translation
       function obj = translate(obj,delta)
           obj.origin = trans3D(obj.origin,delta);
       end
       
       %% Rotation
       function obj = rotX(obj,angleRadians)
           % Rotate around global x-axis
           obj.x_axis = rotx3D(obj.x_axis,angleRadians);
           obj.y_axis = rotx3D(obj.y_axis,angleRadians);
           obj = obj.setZaxis;
       end
       
       function obj = rotY(obj,angleRadians)
           % Rotate around global y-axis
           obj.x_axis = roty3D(obj.x_axis,angleRadians);
           obj.y_axis = roty3D(obj.y_axis,angleRadians);
           obj = obj.setZaxis;
       end
       
       function obj = rotZ(obj,angleRadians)
           % Rotate around global z-axis
           obj.x_axis = rotz3D(obj.x_axis,angleRadians);
           obj.y_axis = rotz3D(obj.y_axis,angleRadians);
           obj = obj.setZaxis;
       end
       
       function obj = rotGRASP(obj,angGRASP)
           % angGRASP = [theta,phi,psi] in radians
           [th,ph,ps] = unpackGRASP(angGRASP);
           % See the GRASP_Technical_Description section 2.1 - all angles in radians
           if numel(angGRASP) > 3
               error('input vector must be of length 3');
           end
           th_hat = obj.x_axis.*cos(th).*cos(ph) + obj.y_axis.*cos(th).*sin(ph) - obj.z_axis.*sin(th);
           ph_hat = -obj.x_axis.*sin(ph) + obj.y_axis.*cos(ph);
           r_hat = obj.x_axis.*sin(th).*cos(ph) + obj.y_axis.*sin(th).*sin(ph) + obj.z_axis.*cos(th);
           obj.x_axis = th_hat.*cos(ph - ps) - ph_hat.*sin(ph - ps);
           obj.y_axis = th_hat.*sin(ph - ps) + ph_hat.*cos(ph - ps);
           obj.z_axis = r_hat;
       end
       
       function obj = rotEuler(obj,angEuler)
           % See the GRASP_Technical_Description section 2.1 - all angles in radians
           if numel(angEuler) > 3
               error('input vector must be of length 3');
           end
           obj = rotGRASP(obj,Euler2GRASP(angEuler));
       end
       
       function Q = dirCosine(coor_new,coor_base)
           % Calculates the 3x3 direction cosine matrix between coordinate
           % systems.  For only one argument, the global coordinate system
           % is assumed as the second system. 
           % coor_new indicates the rotated system and coor_base the base system.
           % So, a point in the rotated system, which was specified in the
           % bases system, is calculated as A_rotate = inv(Q)*A_base
           % See the note: http://homepages.engineering.auckland.ac.nz/~pkel015/SolidMechanicsBooks/Part_III/Chapter_1_Vectors_Tensors/Vectors_Tensors_05_Coordinate_Transformation_Vectors.pdf
           if nargin == 1
               coor_base = coordinateSystem();
           end
           Q = [dot(coor_base.x_axis,coor_new.x_axis), dot(coor_base.x_axis,coor_new.y_axis), dot(coor_base.x_axis,coor_new.z_axis);...
                dot(coor_base.y_axis,coor_new.x_axis), dot(coor_base.y_axis,coor_new.y_axis), dot(coor_base.y_axis,coor_new.z_axis);...
                dot(coor_base.z_axis,coor_new.x_axis), dot(coor_base.z_axis,coor_new.y_axis), dot(coor_base.z_axis,coor_new.z_axis)];
       end
       
       function Uprime = changeBase(U,coor_new,coor_base)
           % Provides the points defined in the [3xN] matrix U, which are
           % defined in the coordinate system coor_base, in the new
           % coordinate system coor_new through translation and rotation.
           % U can also be a coordinateSystem type, the Uprime is returned
           % as the new coordinate system - rotated and translated... 
           % The base class is deleted in this process.
           if nargin == 2
               coor_base = coordinateSystem();
           end
           coorClass = false;
           if isa(U,'coordinateSystem')
               Ncoor = numel(U); % Handle arrays or ccordinate systems
               Utemp = zeros(3,3,Ncoor);
               for ii = 1:Ncoor
                   Utemp(:,:,ii) = [[0;0;0],U(ii).x_axis,U(ii).y_axis]+U(ii).origin;
               end
               U = reshape(Utemp,3,3*Ncoor);
               coorClass = true;
           end
           [N3,~] = size(U);
           if N3 ~= 3
               error('U must have 3 rows indicated [x;y;z] coordinates')
           end
           
           % Move points to origin reference
           Uorigin = bsxfun(@minus,U,coor_base.origin);
           % Rotate the points in the origin reference
           Q = dirCosine(coor_new,coor_base);
           UoriginPrime = Q*Uorigin;
           % Move to new coordinate base
           Uprime = UoriginPrime + coor_new.origin;
           if coorClass
               UprimeTemp = reshape(Uprime,3,3,Ncoor);
               clear Uprime
               Uprime(1,Ncoor) = coordinateSystem(); % Initialise an arry of coordinate systems
               for ii = 1:Ncoor
                   Uprime(ii) = coordinateSystem(UprimeTemp(:,1,ii),UprimeTemp(:,2,ii)-UprimeTemp(:,1,ii),UprimeTemp(:,3,ii)-UprimeTemp(:,1,ii));
               end
           end
       end
       
       function [angGRASP] = getGRASPangBetweenCoors(coor1,coor0)
           % Returns the GRASP angles (in rad) required to rotate from
           % coor0 to coor1.  That is coor1 = coor0.rotGRASP(th,ph,ps)
           % Get the required angles to rotate back to the GC
           % See the GRASP technical description for details on the
           % definitions
           % coor0 corresponds to xyz
           % coor1 corresponds to x1y1z1
           % For only one argument, the global system is assumed for coor0
           
           if nargin == 1
               coor0 = coordinateSystem();
           end
           
           x = coor0.x_axis;
           z = coor0.z_axis;
           x1 = coor1.x_axis;
           z1 = coor1.z_axis;
           
           % Get th - angles between z and z1
           th = angBetweenVectors(z,z1);
           
           % Now get ph
           % Get the vector which is x rotated by ph in the x-y plane 
           Nz = cross(z,z1);
           Nz = Nz./norm(Nz);
           x_ph = cross(Nz,z);
           x_ph = x_ph./norm(x_ph);
           ph = angBetweenVectors(x,x_ph);
           % Sort out the sign of ph - compare to the z-axis direction
           phSign = sign(dot(cross(x,x_ph),z));
           ph = ph*phSign;
           
           % Rotate the system to get to the intermediate coordinate system
           % xyz_prime
           coorPrime = coor0.rotGRASP([th,ph,0]);
           xp = coorPrime.x_axis;
           % Calculate psi
           ps = angBetweenVectors(xp,x1);
           % Sort out the sign of ph - compare to the z1-axis direction
           psSign = sign(dot(cross(xp,x1),z1));
           ps = ps*psSign;
           angGRASP = [th,ph,ps];
       end
       
       % Make a simplified version - this just returns the GRASP angles in
       % the global coordinate system of a given coordinate system
       function angGRASP = getGRASPangles(coor)
           angGRASP = getGRASPangBetweenCoors(coor);
       end
       
       function angEuler = getEulerAngles(coor)
           angEuler = GRASP2Euler(getGRASPangBetweenCoors(coor));
       end
       
       function  B = isequal(coor1,coor2,tol)
            % Test for equality - with a tolerance
            % Does not check the bases...
           if nargin == 2
               tol = norm(coor1.origin).*1e-10;
           end
           BO = all(abs(coor1.origin - coor2.origin) < tol);
           Bx = all(abs(coor1.x_axis - coor2.x_axis) < tol);
           By = all(abs(coor1.y_axis - coor2.y_axis) < tol);
           B = BO && Bx && By;
       end
              
       %% Plotting
       function plot(obj,scale)
           if nargin  == 1
               scale = 1;
           end
           lineWidth = 1;
           textSize = 12;
           x = obj.x_axis.*scale;
           y = obj.y_axis.*scale;
           z = obj.z_axis.*scale;
           xAx = [obj.origin, obj.origin+x];
           yAx = [obj.origin, obj.origin+y];
           zAx = [obj.origin, obj.origin+z];
           plot3(xAx(1,:),xAx(2,:),xAx(3,:),'r','LineWidth',lineWidth), hold on
           plot3(yAx(1,:),yAx(2,:),yAx(3,:),'g','LineWidth',lineWidth)
           plot3(zAx(1,:),zAx(2,:),zAx(3,:),'b','LineWidth',lineWidth)
           text(xAx(1,2),xAx(2,2),xAx(3,2),'x','FontSize',textSize,'color','r')
           text(yAx(1,2),yAx(2,2),yAx(3,2),'y','FontSize',textSize,'color','g')
           text(zAx(1,2),zAx(2,2),zAx(3,2),'z','FontSize',textSize,'color','b')
           grid on
           axis equal
           xlabel('x-axis (m)')
           ylabel('y-axis (m)')
           zlabel('z-axis (m)')
           view([140,40])
       end
       
   end
   
   methods (Access = private)
       function obj = setZaxis(obj)
           obj.z_axis = getZaxis(obj);
       end
       function z = getZaxis(obj)
           z = cross(obj.x_axis,obj.y_axis);
       end
       function obj = normAxis(obj)
           % Make sure we have unit vectors
           obj.x_axis = obj.x_axis/norm(obj.x_axis);
           obj.y_axis = obj.y_axis/norm(obj.y_axis);
           obj.z_axis = obj.z_axis/norm(obj.z_axis);
       end
   end
end