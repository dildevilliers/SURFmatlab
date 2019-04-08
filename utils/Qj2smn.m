% function [Qsmn] = Qj2smn(Qj)
% Converts j-indexed Q-coefficient struct (as output by FF2SWE.m) to its smn-indexed
% counterpart.
% Inputs:
%   --Qj: Struct as output by FF2SWE.m, containing SWE Q-coefficients using
%   single j-indexing (as in Hansen)
% Returns:
%   --Qsmn: Struct as output by FF2SWE.m, containing SWE Q-coefficients using
%   default smn-indexing (as in Hansen)
%NB: See smn2j.m and j2smn.m for quick conversions between the two indexing
%schemes.

% Created: 2017-10-24
% Brandt Klopper


function Qsmn = Qj2smn(Qj)

Qsmn = rmfield(Qj,{'Q'});

for xx = 1:size(Qj,1)
    for ff = 1:size(Qj,2)
    
        %Build Qsm0n
        for ss = 1:2
            for nn = 1:Qj(xx,ff).NMAX
                J = smn2j([ss,0,nn]);
                Qsmn(xx,ff).Qsm0n(ss,1,nn) = Qj(xx,ff).Q(J);
            end
        end

        [Qsmn(xx,ff).Qsmmn,Qsmn(xx,ff).Qsmpn] = deal(zeros(2,Qj(xx,ff).MMAX,Qj(xx,ff).NMAX));

        %Build Qsmmn and Qsmpn
        for ss = 1:2
            for mm = 1:Qj(xx,ff).MMAX
                for nn = 1:Qj(xx,ff).NMAX
                    if mm <= nn
                        Jm = smn2j([ss,-mm,nn]);
                        Jp = smn2j([ss,+mm,nn]);
                        Qsmn(xx,ff).Qsmmn(ss,mm,nn) = Qj(xx,ff).Q(Jm);
                        Qsmn(xx,ff).Qsmpn(ss,mm,nn) = Qj(xx,ff).Q(Jp);
                    end
                end
            end
        end
    
    end
end

end