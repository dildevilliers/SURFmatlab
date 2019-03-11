function smn = j2smn(J)

if mod(J,2) ~= 0
    s = 1;
else
    s = 2;
end

n = floor(sqrt((J-s)/2 + 1)); 

m = (J-s)/2 + 1 - n*(n+1);

smn = [s;m;n];

end