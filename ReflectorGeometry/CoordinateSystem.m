classdef CoordinateSystem
    properties
      origin(1,1) Pnt3D = Pnt3D(0,0,0) 
      base = [] % Can be another coordinate system object, or, if empty, is assumed to be the global coordinate system 
    end
   
    properties (SetAccess = private)
        x_axis = [1;0;0] % x-axis direction
        y_axis = [0;1;0] % y-axis direction
    end
    
    properties (Dependent = true)
        z_axis
    end
    
    methods
       function obj = CoordinateSystem(origin,x_axis,y_axis,base)
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
               if ~(isempty(base.base))
                   error('base coordinate system must be based on the global coordinates - use getInGlobal to get base coordinates in global coordinates')
               else
                   obj.base = base;
               end
           end
           % Check for valid inputs
           nX = norm(obj.x_axis);
           nY = norm(obj.y_axis);
           z = cross(obj.x_axis,obj.y_axis)/(nX*nY);
           if abs(norm(z)-1) > 1e-10
               error('x_axis and y_axis must be orthogonal');
           end
           % Set the z-axis direction
           % Normalise the unit vectors
           obj = obj.normAxis;
       end
       
       %% Setters
       function z_axis = get.z_axis(obj)
            z_axis = cross(obj.x_axis,obj.y_axis);
       end
       
       function obj = set2Base(obj)
           if isa(obj.base,'CoordinateSystem')
               obj = obj.base;
           end
       end
       
       %% Translation
       function obj = translate(obj,delta)
           obj.origin = addVect(obj.origin,delta);
%            obj.origin = translate(obj.origin,delta);
       end
       
       %% Rotation
       function obj = rotX(obj,angleRadians)
           % Rotate around global x-axis
           obj.x_axis = rotx(obj.x_axis,angleRadians);
           obj.y_axis = rotx(obj.y_axis,angleRadians);
       end
       
       function obj = rotY(obj,angleRadians)
           % Rotate around global y-axis
           obj.x_axis = roty(obj.x_axis,angleRadians);
           obj.y_axis = roty(obj.y_axis,angleRadians);
       end
       
       function obj = rotZ(obj,angleRadians)
           % Rotate around global z-axis
           obj.x_axis = rotz(obj.x_axis,angleRadians);
           obj.y_axis = rotz(obj.y_axis,angleRadians);
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
           % is assumed as the second (base) system. 
           % coor_new indicates the rotated system and coor_base the base system.
           % So, a point in the rotated system, which was specified in the
           % base system, is calculated as A_rotate = inv(Q)*A_base
           % See the note: http://homepages.engineering.auckland.ac.nz/~pkel015/SolidMechanicsBooks/Part_III/Chapter_1_Vectors_Tensors/Vectors_Tensors_05_Coordinate_Transformation_Vectors.pdf
           if nargin == 1
               coor_base = CoordinateSystem();
           end
           Q = [dot(coor_base.x_axis,coor_new.x_axis), dot(coor_base.x_axis,coor_new.y_axis), dot(coor_base.x_axis,coor_new.z_axis);...
                dot(coor_base.y_axis,coor_new.x_axis), dot(coor_base.y_axis,coor_new.y_axis), dot(coor_base.y_axis,coor_new.z_axis);...
                dot(coor_base.z_axis,coor_new.x_axis), dot(coor_base.z_axis,coor_new.y_axis), dot(coor_base.z_axis,coor_new.z_axis)];
       end
       
       %% Change of basis
       function coorInGlobal = getInGlobal(obj)
           % Returns the coordinate system, defined in obj.base, in the global coordinate system
           % Assume obj.base is always defined in the global coordinate system
           
           if ~isempty(obj.base)
               oldBase = obj.base;
               newBase = CoordinateSystem();
               % Since the coordinate system is always defined in the global
               % coordinate system, with a given base (which corresponds to the
               % local global coordinate system), get the angle between the
               % new base coordinate system and the object
               graspAng = getGRASPangles(obj);
               % Now rotate the base by this amount
               baseRotated = oldBase.rotGRASP(graspAng);
               % And shift the origin
               Ups = changeBase(obj.origin,newBase,oldBase);
               coorInGlobal = CoordinateSystem(Ups,baseRotated.x_axis,baseRotated.y_axis);
               coorInGlobal.base = [];
           else
               coorInGlobal = obj;
           end
       end
       
       function coorOut = redefineToOtherBase(obj,newBase)
           % Redefines the coordinate system in obj to be relative to the
           % new base coordinate system newBase
           
           % First get both coordinate systems in the global base
           coorIn = obj.getInGlobal;
           baseIn = newBase.getInGlobal;
           % Calculate the translation
           diffGlobal = coorIn.origin.pointMatrix - baseIn.origin.pointMatrix;
           Ox = dot(diffGlobal,baseIn.x_axis);
           Oy = dot(diffGlobal,baseIn.y_axis);
           Oz = dot(diffGlobal,baseIn.z_axis);
           % Calculate the rotation
           graspAng = getGRASPangBetweenCoors(coorIn,baseIn);
           % Build the coordinate system
           coorOut = CoordinateSystem(Pnt3D(Ox,Oy,Oz));
           coorOut = coorOut.rotGRASP(graspAng);
           coorOut.base = newBase;
       end
       
       function [angGRASP] = getGRASPangBetweenCoors(coor1,coor0)
           % Returns the GRASP angles (in rad) required to rotate from
           % coor0 to coor1.  That is coor1 = coor0.rotGRASP([th,ph,ps])
           % See the GRASP technical description for details on the
           % definitions
           % coor0 corresponds to xyz
           % coor1 corresponds to x1y1z1
           % For only one argument, the global system is assumed for coor0
           
           if nargin == 1
               coor0 = CoordinateSystem();
           end
           
           x = coor0.x_axis;
           z = coor0.z_axis;
           x1 = coor1.x_axis;
           z1 = coor1.z_axis;
           
           % Get th - angles between z and z1
           th = angBetweenVectors(z,z1);
           
           % Now get ph
           % Get the vector which is x rotated by ph in the x-y plane 
           if ~isequal(abs(z),abs(z1))
               Nz = cross(z,z1);
               Nz = Nz./norm(Nz);
               x_ph = cross(Nz,z);
               x_ph = x_ph./norm(x_ph);
           else
               x_ph = x1;
           end
           ph = angBetweenVectors(x,x_ph);
           % Sort out the sign of ph - compare to the z-axis direction
           phSign = sign(dot(cross(x,x_ph),z));
           if phSign ~= 0
               ph = ph*phSign;
           end
           
           % Rotate the system to get to the intermediate coordinate system
           % xyz_prime
           coorPrime = coor0.rotGRASP([th,ph,0]);
           xp = coorPrime.x_axis;
           % Calculate psi
           ps = angBetweenVectors(xp,x1);
           % Sort out the sign of ps - compare to the z1-axis direction
           psSign = sign(dot(cross(xp,x1),z1));
           if psSign ~= 0
               ps = ps*psSign;
           end
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
       
       %% Testers
       function  B = isequal(coor1,coor2,tol)
            % Test for equality - with a tolerance
            % Does not check the bases...
           if nargin == 2
               tol = coor1.origin.r.*1e-10;
           end
           BO = isequal(coor1.origin,coor2.origin);
           Bx = all(abs(coor1.x_axis - coor2.x_axis) < tol);
           By = all(abs(coor1.y_axis - coor2.y_axis) < tol);
           B = BO && Bx && By;
       end
              
       %% Plotting
       function plot(obj,scale)
           % plots the coordinate system in obj in the global coordinates
           if nargin  == 1
               scale = 1;
           end
           coorGlobalBase = obj.getInGlobal;
           plotLocal(coorGlobalBase,scale);
       end
       
       function plotLocal(obj,scale)
           % plots the coordinate system in obj in the local base
           % coordinates
           if nargin  == 1
               scale = 1;
           end
           lineWidth = 1;
           textSize = 12;
           x = obj.x_axis.*scale;
           y = obj.y_axis.*scale;
           z = obj.z_axis.*scale;
           O = [obj.origin.x;obj.origin.y;obj.origin.z];
           xAx = [O, O+x];
           yAx = [O, O+y];
           zAx = [O, O+z];
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
       function obj = normAxis(obj)
           % Make sure we have unit vectors
           obj.x_axis = obj.x_axis/norm(obj.x_axis);
           obj.y_axis = obj.y_axis/norm(obj.y_axis);
       end
   end
end
