classdef (Abstract) FarFieldExpansion
    properties(Abstract, SetAccess = private)
        nBasis
        basis
        nCoeffs
        coeffs
    end
    
    methods
        function plotCoeffs(obj,varargin)
            
            p = inputParser;
            addParameter(p,'plotAxis','coeffs');
            addParameter(p,'plotType','mag');
            addParameter(p,'freq',obj.freq);
            addParameter(p,'coeff',1:obj.nCoeffs);
            addParameter(p,'par',1:length(obj.coeffs));
            
            parse(p, varargin{:})
            plotAxis = p.Results.plotAxis;
            plotType = p.Results.plotType;
            freq = p.Results.freq;
            coeff= p.Results.coeff;
            par = p.Results.par;
            
            ifreq = find(ismember(obj.freq,freq));
            nn = 1;
            
            switch plotAxis
                case 'coeffs'
                    x = 1:obj.nCoeffs;
                    for hh = 1:length(par)
                        for ii = 1:length(freq)
                            for jj = 1:1%length(coeff)
                                y(:,nn) = obj.coeffs{par(hh)}(:,ifreq(ii));
                                plotleg{nn} = ['par.',num2str(par(hh)),'-',num2str(freq(ii)),' Hz'];
                                nn = nn + 1;
                            end
                        end
                    end
                case 'freq'
                    x = obj.freq;
                    for hh = 1:length(par)
                        for ii = 1:1%length(freq)
                            for jj = 1:length(coeff)
                                y(:,nn) = obj.coeffs{par(hh)}(coeff(jj),:).';
                                plotleg{nn} = ['par.',num2str(par(hh)),'-C',num2str(coeff(jj))];
                                nn = nn + 1;
                            end
                        end
                    end
                otherwise
                    error('plotAxis string unrecognised')
            end
            
            figure
            switch plotType
                case 'mag'
                    semilogy(x,abs(y));
                    grid on
                case 'pha'
                    semilogy(x,ang(y));
                    grid on
                case 're'
                    plot(x,real(y));
                    grid on
                case 'im'
                    plot(x,imag(y));
                    grid on
                case 'magpha'
                case 'reim'
                otherwise
            end
            legend(plotleg);
            
        end
        
        
        
        function interpCoeffs(obj,coeffs,absc,plotType)
            %Write me!
        end
    end
    
    methods (Abstract, Static)
        coeffs2FarField
        farField2Coeffs
        %getBasisPower
    end
    
end