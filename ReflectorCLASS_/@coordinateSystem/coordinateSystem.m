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
       
       function obj = rotGRASP(obj,theta,phi,psi)
           % See the GRASP_Technical_Description section 2.1 - all angles in radians
           if numel(theta) > 1 || numel(phi) > 1 || numel(psi) > 1
               error('Only scalar input angles accepted');
           end
          th_hat = obj.x_axis.*cos(theta).*cos(phi) + obj.y_axis.*cos(theta).*sin(phi) - obj.z_axis.*sin(theta);
          ph_hat = -obj.x_axis.*sin(phi) + obj.y_axis.*cos(phi);
          r_hat = obj.x_axis.*sin(theta).*cos(phi) + obj.y_axis.*sin(theta).*sin(phi) + obj.z_axis.*cos(theta);
          obj.x_axis = th_hat.*cos(phi - psi) - ph_hat.*sin(phi - psi);
          obj.y_axis = th_hat.*sin(phi - psi) + ph_hat.*cos(phi - psi);
          obj.z_axis = r_hat;
       end
       
       function obj = rotEuler(obj,alpha,beta,gamma)
           % See the GRASP_Technical_Description section 2.1 - all angles in radians
           if numel(alpha) > 1 || numel(beta) > 1 || numel(gamma) > 1
               error('Only scalar input angles accepted');
           end
           theta = beta;
           phi = alpha;
           psi = alpha + gamma;
           obj = rotGRASP(obj,theta,phi,psi);
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
           % defined in the coordiante system coor_base, in the new
           % coordinate system coor_new through translation and rotation.
           if nargin == 2
               coor_base = coordinateSystem();
           end
           [N3,~] = size(U);
           if N3 ~= 3
               error('U must have 3 rows indictaed [x;y;z] coordinates')
           end
           % Move points to origin reference
           Uorigin = U - coor_base.origin;
           % Raotate the points in the origin reference
           Q = dirCosine(coor_new,coor_base);
           UoriginPrime = Q*Uorigin;
           % Move to new coordinate base
           Uprime = UoriginPrime + coor_new.origin;
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