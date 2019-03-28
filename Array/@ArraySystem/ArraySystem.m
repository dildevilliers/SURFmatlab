classdef ArraySystem

    properties
        % Elements
        antPos(1,1) pnt3D = pnt3D([0,1],0,0)   % Antenna positions in 3D points - internal vector
        channelPhasors(1,:) double {mustBeFinite} = 1 % Vector of complex channel phasors for scaling/calibration
        % Receiver
        noisePower(1,1) double {mustBeReal, mustBeFinite} = -110 % at the input and in dBm
        LNAGain(1,1) double {mustBeReal, mustBeFinite} = 20    % LNA gain in dB
        IFGain(1,1) double {mustBeReal, mustBeFinite} = 130    % IF gain in dB
        freqLO(1,1) double {mustBeReal, mustBeFinite} = 1.571324e9        % LO frequency in Hz
        % Timing
        freqSamp(1,1) double {mustBeReal, mustBeFinite} = 16.368e6   % Sample rate in Hz
        Nt(1,1) double {mustBeReal, mustBeFinite} = 1600;            % Number of time samples
        % ADC
        Nbits(1,1) double {mustBeInteger, mustBeFinite} = 1 
        maxV(1,1) double {mustBeReal, mustBeFinite} = 1
        minV(1,1) double {mustBeReal, mustBeFinite} = -1
        ADCtype(1,:) char {mustBeMember(ADCtype,{'binary','sign_mag','2scomp'})} = 'binary' 
    end
    
    properties (SetAccess = private)
        elements
        receiver
        adc
        delT
        t
    end
    
    methods
         % Constructor
         function obj = ArraySystem(antPos,channelPhasors,noisePower,LNAGain,IFGain,freqLO,freqSamp,Nt,Nbits,maxV,minV,ADCtype)
             if nargin >= 1, obj.antPos = antPos; end
             if nargin >= 2, obj.channelPhasors = channelPhasors; end
             if nargin >= 3, obj.noisePower = noisePower; end
             if nargin >= 4, obj.LNAGain = LNAGain; end
             if nargin >= 5, obj.IFGain = IFGain; end
             if nargin >= 6, obj.freqLO = freqLO; end
             if nargin >= 7, obj.freqSamp = freqSamp; end
             if nargin >= 8, obj.Nt = Nt; end
             if nargin >= 9, obj.Nbits = Nbits; end
             if nargin >=10, obj.maxV = maxV; end
             if nargin >=11, obj.minV = minV; end
             if nargin >=12, obj.ADCtype = ADCtype; end
             
             obj.elements = ArrayElements(obj.antPos,obj.channelPhasors);
             obj.receiver = ArrayReceiver(obj.noisePower,obj.LNAGain,obj.IFGain,obj.freqLO);
             obj.adc = ArrayADC(obj.Nbits,obj.maxV,obj.minV,obj.ADCtype);
             obj.delT = 1/obj.freqSamp;
             t0 = 0;   % Hardcode for now - no reason to change I think...
             obj.t = t0:obj.delT:(t0+obj.delT*(obj.Nt-1));
         end
         
         function x = getPortSignal(obj,s,Qtype)
             % Returns the digitised signals at the antenna ports as a
             % matrix of size [Nant, Nt].  s is an array of PlaneWaveSignal
             % objects
             % Qtype describes how to handle the Q channel:
             %  0: use the actual Q
             %  [...-2,-1,1,2,...]: Integer shift by this number
             %  +-inf: Hilbert transform
             
             if nargin < 3
                 Qtype = 0;
             end
             assert(isscalar(Qtype),'Error: Qtype must be scalar');
             
             % Get signals after elements
             portSigMat = obj.elements.portSignals(s,obj.t);
             % Put them through receiver
             sn = obj.receiver.sigRec(portSigMat,obj.t);
             % Digitise
             xi = obj.adc.ADC(real(sn));
             if Qtype == 0
                 xq = obj.adc.ADC(imag(sn));
             elseif isinteger(Qtype)
                 xq = circshift(xi,Qtype,2);
             elseif isinf(Qtype)
                 xq = hilbert(xi);
             else
                 error(['Unknown Qtype'])
             end
             x = xi + 1i.*xq;
         end
         
         function plotPortSignal(obj,s,portNumber)
             if nargin < 3
                 portNumber = 1;
             end
             x = obj.getPortSignal(s);
             xP = real(x(portNumber,:));
             plot(obj.t,xP), grid on, hold on
             xlabel('t (s)')
         end
    end
    
end
