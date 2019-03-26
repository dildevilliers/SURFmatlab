classdef ArrayElements
    %ARRAYELEMENTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        antPos(1,1) pnt3D = pnt3D([0,1],0,0)   % Antenna positions in 3D points - internal vector
        channelErrors(1,:) double {mustBeFinite} = 1 % Vector of complex channel errors for calibration testing
    end
    
    properties (SetAccess = private)
       N_elements     % Number of elements
    end
    
    methods
        function obj = ArrayElements(antPos,channelErrors)
            if nargin >= 1
                obj.antPos = antPos;
            end
            if nargin >= 2
                obj.channelErrors = channelErrors;
            end
            obj.N_elements = length(obj.antPos.x);
            if length(obj.channelErrors) == 1
                obj.channelErrors = repmat(obj.channelErrors,1,obj.N_elements);
            else 
                assert(length(obj.channelErrors) == obj.N_elements,'Length of channelErrors should match the number of elements');
            end
        end
        
        function portSigMat = portSignals(obj,S,t)
            % Takes an array of PlaneWaveSignal objects (in S), adds them up, and
            % calculates the complex analog signals at the array antenna ports  at the time samples t 
            if nargin == 1
                S = PlaneWaveSignal;
                t = linspace(0,10,1001);
            end
            assert(min(size(t)) == 1,'Error, t must be a vector, not a matrix');
            assert(isreal(t),'Error, t must be real');
            t = t(:);
            Nsig = length(S);
            % Build the signal matrix
            si = zeros(length(t),Nsig);
            k = zeros(3,Nsig);
            for ss = 1:Nsig
                si(:,ss) = S(ss).getSignal(t);
                k(:,ss) = S(ss).k;
            end
            % Get the element phase shifts
            r = obj.antPos.pointMatrix.';  % [Nant x 3]
            A = exp(-1i*r*k);              % [Nant x Nsig]
            
            % Add the systematic errors in the channels
            AC = bsxfun(@times,A,obj.channelErrors(:));
            
            % Multiply the signals through the channels
            portSigMat = AC*si.';
        end
    end
end

