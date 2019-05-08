classdef ellipticalRim
    % Info found:
    % https://math.stackexchange.com/questions/426150/what-is-the-general-equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate
    % https://math.stackexchange.com/questions/2645689/what-is-the-parametric-equation-of-a-rotated-ellipse-given-the-angle-of-rotatio?noredirect=1&lq=1
   properties
      centre(2,1) double {mustBeReal, mustBeFinite} = [0;0] % [x,y] in (m) of rim centre
      halfAxis(2,1) double {mustBeReal, mustBeFinite} = [1;1] % [a,b] in (m) - half axis along x = a, and along y = b 
      rotation(1,1) double {mustBeReal, mustBeFinite} = 0 % rotation angle in (rad) from x-axis
   end
   methods
       function obj = ellipticalRim(cen,hA,rot)
           if nargin == 0
           elseif nargin == 1
               obj.centre = cen;
           elseif nargin == 2
               obj.centre = cen;
               obj.halfAxis = hA;
           else
               obj.centre = cen;
               obj.halfAxis = hA;
               obj.rotation = rot;
           end
       end
       
       function M = isInRim(obj,x,y)
           % Test if the point [x,y] is in the rim - return boolean values
           [Rx,Ry,Cx,Cy,alpha] = unpackEll(obj);
           M = ((x-Cx).*cos(alpha) + (y-Cy).*sin(alpha)).^2./Rx.^2 + ((x-Cx).*sin(alpha) - (y-Cy).*cos(alpha)).^2./Ry.^2;
           M(M<=1) = 1;
           M(M>1) = 0;
       end
       
       function V = cartRim(obj,N)
           if nargin == 1
               N = 101;
           end
           % Return cartesian vectors, of size [2xN], along the rim [x,y].'
           [Rx,Ry,Cx,Cy,alpha] = unpackEll(obj);
           ph = linspace(0,2*pi,N);
           x = Cx + Rx.*cos(ph).*cos(alpha) - Ry.*sin(ph).*sin(alpha);
           y = Cy + Rx.*cos(ph).*sin(alpha) + Ry.*sin(ph).*cos(alpha);
%            V = [x;y];
           V = Pnt3D(x,y,0);
       end
       
       function V = polarRim(obj,N)
           if nargin == 1
               N = 101;
           end
           % Return polar vectors, of size [2xN], along the rim [ph,rho].'
           Vc = cartRim(obj,N);
           [ph,rho] = cart2pol(Vc(1,:),Vc(2,:));
           V = [ph;rho];
       end
       
       function plot(obj)
           V = cartRim(obj,100);
%            plot(V(1,:),V(2,:),'k')
           plot(V.x,V.y,'k')
           grid on
           xlabel('x-axis (m)')
           ylabel('y-axis (m)')
       end
       
   end
   
   methods (Access = private)
       function [Rx,Ry,Cx,Cy,alpha] = unpackEll(obj)
           Rx = obj.halfAxis(1);
           Ry = obj.halfAxis(2);
           Cx = obj.centre(1);
           Cy = obj.centre(2);
           alpha = obj.rotation;
       end
   end
end