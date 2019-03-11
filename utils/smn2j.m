function J = smn2j(smn)

[s,m,n] = deal(smn(1),smn(2),smn(3));

J = 2*( n*(n+1) + m - 1) + s;

end