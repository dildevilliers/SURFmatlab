% function [Qj] = Qsmn2j(Qsmn)
% Converts smn-indexed Q-coefficient struct (as output by FF2SWE.m and readSPH.m) to its j-indexed
% counterpart.
% Inputs:
%   --Qsmn: Struct as output by FF2SWE.m, containing SWE Q-coefficients using
%   default smn-indexing (as in Hansen)
% Returns:
%   --Qj: Struct as output by FF2SWE.m, containing SWE Q-coefficients using
%   single j-indexing (as in Hansen)
%NB: See smn2j.m and j2smn.m for quick conversions between the two indexing
%schemes.

% Created: 2017-10-24
% Brandt Klopper

function Qj = Qsmn2j(Qsmn)

Qj = rmfield(Qsmn,{'Qsm0n','Qsmmn','Qsmpn'});

for xx = 1:size(Qsmn,1)
    for ff = 1:size(Qsmn,2)
    
        Qj(xx,ff).Q = zeros(2*Qsmn(xx,ff).NMAX*(Qsmn(xx,ff).NMAX+2),1);

        for ss = 1:2
            for mm = -Qsmn(xx,ff).MMAX:Qsmn(xx,ff).MMAX
                for nn = 1:Qsmn(xx,ff).NMAX
                    if abs(mm) <= nn
                        J = smn2j([ss mm nn]);
                        if mm == 0
                            Qj(xx,ff).Q(J) = Qsmn(xx,ff).Qsm0n(ss,1,nn);
                        elseif mm < 0
                            Qj(xx,ff).Q(J) = Qsmn(xx,ff).Qsmmn(ss,-mm,nn);
                        else
                            Qj(xx,ff).Q(J) = Qsmn(xx,ff).Qsmpn(ss,mm,nn);
                        end
                    end
                end
            end
        end

    end
end

end